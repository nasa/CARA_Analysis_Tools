function [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples,Frequency,SwerlingType,Ln_MassCalFac_Mean,Ln_MassCalFac_Sigma)
%
% EstimateObjectMassQuantiles - Estimates a Secondary Objects Mass on a quantile level from an
%                               input RCS, Ballistic Coefficient, Ballistic
%                               Coefficient Variance, Drag Coefficient, and
%                               Drag Coefficient Variance (or limits), These are
%                               determined by generating a large number of
%                               samples in a Monte Carlo Fashion and
%                               reports values out on a quantile basis
%                               (The code can also work if instead of Drag
%                               coefficients and ballistic coefficients,
%                               Reflectivity coefficient and SRP
%                               coefficients are provided)
%
% Syntax:   
%   [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar)
%   [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar)
%   [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector)
%   [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples)
%   [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples,Frequency)
%   [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples,Frequency,SwerlingType)
%   [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples,Frequency,SwerlingType,Ln_MassCalFac_Mean,Ln_MassCalFac_Sigma)
%
% Inputs:
%
%   RCS             - 1X1 Radar Cross Section of Secondary Object (m^2)
%   B               - 1X1 Ballistic Coefficient of Secondary Object (m^2/kg)
%   BVar            - 1X1 Ballistic Coefficient Variance of Secondary
%                     Object (m^4/kg^2)
%   Cd              - 1X1 Drag Coefficient Estimate of Secondary Object (dimensionless)
%                     (optional, default: 2.5)
%   CdVar           - 1X1 or 1X2 Drag Coefficient Variance of Secondary
%                     Object if only 1 element is provided. If 2 elements are 
%                     provided, Cd vector will be sampled uniformly instead
%                     of normally, where the lower bound and maximum bounds
%                     are Cd-CdVar(1) and Cd+CdVar(2) respectively
%                     (optional, default: [0.4 0.4])
%   QuantileVector  - 1XN Vector of desired mass quantile estimates
%                     (optional, defaults to 0.999 (99.9%))
%   NumOfSamples    - INTEGER number of samples to generate (optional,
%                     Default = 1E5, or a minimum number to accurately describe highest input quantile)
%   Frequency       - 1X1 Frequency of the sensor at which RCS have been
%                     obtained (MHz) (optional, Default = 2000 MHz)
%   SwerlingType    - Text input of Swerling distribution type (optional, Default = 'III')
%                     Allowable Inputs:
%                       * 'I'
%                       * 'II'
%                       * 'III'
%                       * 'IV'
%   Ln_MassCalFac_Mean - 1X1 The mean log-normal mass calibration factor to
%                        scale the calculated masses (optional, Default = 0)
%   Ln_MassCalFac_Sigma- 1X1 The std log-normal mass calibration factor to
%                        scale the calculated masses (optional, Default = 0)
%
% Outputs:
%   massVec         - NumOfSamplesX1 array of the secondary object mass
%                     estimates for each individual sample
%   QuantileArray   - NX2 Structure array of Mass Quantile Estimates
%   AreaVec         - NumOfSamplesX1 Estimated frontal area vector assuming 
%                     a spherical cross section
%
%
% Other m-files required: 	RCSDistribution.m
%                           EstimateMassFromRCS.m
% Subfunctions: None
% MAT-files required:       None
%
% See also: Lechtenberg, T., "An Operational Algorithm for Evaluating
%           Satellite Collision Consequence," AAS Astrodynamics Specialist
%           Conference, 2019, AAS 19-669
%
% May 2019; Last revision: 04-Jun-2026
%
% ----------------- BEGIN CODE -----------------
    
    %% Set up Defaults
    % Default Drag coefficient
    if nargin < 4 || isempty(Cd)
        Cd = 2.5;
    end
    if nargin < 5 || isempty(CdVar)
        CdVar = [0.4 0.4];
    end

    % Default Mass Quantile
    if nargin < 6 || isempty(QuantileVector)
        QuantileVector = 0.999;
    end
    
    % Default number of samples
    if nargin < 7 || isempty(NumOfSamples)
        NumOfSamples = round(max(1e5,10/(1-max(QuantileVector))));
    end

    % Default Frequency
    if nargin < 8 || isempty(Frequency)
        Frequency = 2000;
    end

    % Set Default RCS Sampling inputs
    if nargin < 9 || isempty(SwerlingType)
        SwerlingType            = 'III'; % Identification of Swerling Gamma Distribution Type
    end

    % Default log-normal mass calibration factor
    if nargin < 11 || isempty(Ln_MassCalFac_Mean) || isempty(Ln_MassCalFac_Sigma)
        Ln_MassCalFac_Mean = 0;
        Ln_MassCalFac_Sigma = 0;
    end
    
    % Drag Coefficient Distribution
    if isscalar(CdVar)
        % normal sampling
        CdSamps=abs(Cd+sqrt(CdVar)*randn(NumOfSamples,1));
    elseif numel(CdVar)==2
        % uniform sampling
        CdLow = Cd-CdVar(1);
        CdSamps = CdLow + (CdVar(1) + CdVar(2))*rand(NumOfSamples,1);
    end
    
    % Sample RCS distribution from Swerling III distribution
    [RCSSamps] = RCSDistribution(RCS,NumOfSamples,SwerlingType);
    
    % ballistic coefficient distribution--normal sampling
    BSamps=abs(B+sqrt(BVar)*randn(NumOfSamples,1));
    
    % Generate a vector of possible Secondary Object Masses
    [massVec,~,AreaVec] = EstimateMassFromRCS(RCSSamps,CdSamps,BSamps,SwerlingType,...
        Frequency,NumOfSamples,Ln_MassCalFac_Mean,Ln_MassCalFac_Sigma);
    
    % Sort Mass Estimation Vector
    SortedMassVec = sort(massVec);

    % Get Quantile Data
    for i = 1:length(QuantileVector)
        QuantileArray(i).EstimationQuantile = QuantileVector(i);
        QuantileArray(i).MassEstimate       = SortedMassVec(round(QuantileVector(i)*length(massVec)));
    end

    

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 05-20-2019 | Initial Development
% S. Es haghi    | 06-04-2026 | Adding more flexibility to the function by
%                               taking additional inputs and correcting
%                               header
%
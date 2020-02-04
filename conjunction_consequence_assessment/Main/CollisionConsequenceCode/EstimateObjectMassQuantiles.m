function [massVec,QuantileArray,AreaVec] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples)
%
% EstimateObjectMassQuantiles - Estimates a Secondary Objects Mass on a quantile level from an
%                               input RCS, Ballistic Coefficient, Ballistic
%                               Coefficient Variance, Drag Coefficient, and
%                               Drag Coefficient Variance, These are
%                               determined by generating a large number of
%                               samples in a Monte Carlo Fashion and
%                               reports values out on a quantile basis
%
% Syntax:   [Catastrophic,NumOfPieces] = CollisionConsequenceNumPieces(PrimaryMass,VRel,RCS,B,BVar,NumOfSamples)
%
% Inputs:

%   RCS             - 1X1 Radar Cross Section of Secondary Object (m^2)
%   B               - 1X1 Ballistic Coefficient of Secondary Object (m^2/kg)
%   BVar            - 1X1 Ballistic Coefficient Variance of Secondary
%                     Object (m^4/kg^2)
%   B               - 1X1 Drag Coefficient Estimate of Secondary Object (dimensionless)
%   BVar            - 1X1 Drag Coefficient Variance of Secondary
%                     Object (dimensionless)
%   QuantileVector  - 1XN Vector of desired mass quantile estimates
%                     (optional, defaults to 0.999 (99.9%))
%   NumOfSamples    - INTEGER number of samples to generate (optional,
%                     Default = 10000, or a minimum number to accurately describe highest input quantile)
%
% Outputs:
%   massVec         - [NumOfSamplesX1] array of the secondary object mass
%                     estimates for each individual sample
%   QuantileArray   - [NX2] Structure array of Mass Quantile Estimates
%   AreaVec         - 
%
% Example/Validation Cases:
%
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
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
% May 2019; Last revision: 20-May-2019
%
% ----------------- BEGIN CODE -----------------
    
    %% Set up Defaults
    % Default Mass Quantile
    if nargin < 6 || isempty(QuantileVector)
        QuantileVector = 0.999;
    end
    
    % Default number of samples
    if nargin < 7 || isempty(NumOfSamples)
        NumOfSamples = round(max(10000,10/(1-max(QuantileVector))));
    end
    
    % Set Default RCS Sampling inputs
    SwerlingType            = 'III'; % Identification of Swerling Gamma Distribution Type
    
    % Drag Coefficient Distribution--normal sampling  
    CdSamps=abs(Cd+sqrt(CdVar)*randn(NumOfSamples,1)); 
    
    % Sample RCS distribution from Swerling III distribution
    [RCSSamps] = RCSDistribution(RCS,NumOfSamples,SwerlingType);
    
    % ballistic coefficient distribution--normal sampling
    BSamps=abs(B+sqrt(BVar)*randn(NumOfSamples,1));
    
    % Generate a vector of possible Secondary Object Masses
    [massVec,~,AreaVec] = EstimateMassFromRCS(RCSSamps,CdSamps,BSamps,SwerlingType);
    
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
%
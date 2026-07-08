function [EstimatedMass,RCS,AreaVec] = ...
        EstimateMassFromRCS(RCS,CdVec,BVec,SwerlingType,Frequency,NumOfSamples,Ln_MassCalFac_Mean,Ln_MassCalFac_Sigma)
%
% EstimateMassFromRCS - Estimates Mass from an input satellite's RCS , 
%                       Ballistic Coefficient and Drag Coefficients
%                       (The code can also work if instead of Drag
%                       coefficients and ballistic coefficients,
%                       Reflectivity coefficient and SRP
%                       coefficients are provided)
%
% Syntax:   
% 
% [EstimatedMass,RCS,AreaVec] = EstimateMassFromRCS(RCS,CdVec,BVec)
% [EstimatedMass,RCS,AreaVec] = EstimateMassFromRCS(RCS,CdVec,BVec,SwerlingType)
% [EstimatedMass,RCS,AreaVec] = EstimateMassFromRCS(RCS,CdVec,BVec,SwerlingType,Frequency)
% [EstimatedMass,RCS,AreaVec] = EstimateMassFromRCS(RCS,CdVec,BVec,SwerlingType,Frequency,NumOfSamples)
% [EstimatedMass,RCS,AreaVec] = EstimateMassFromRCS(RCS,CdVec,BVec,SwerlingType,Frequency,NumOfSamples,Ln_MassCalFac_Mean,Ln_MassCalFac_Sigma)
%
% Inputs:
%   RCS             - NX1 or 1X1 Radar Cross Section of Object (m^2)
%   CdVec           - NX1 or 1X1 or 1X2 Coefficient of Drag of the object (dimensionless)
%                     (if a 2 element vector is provided, the first element
%                     will be taken as coefficient mean and the second
%                     element will be taken as coefficient variance and a
%                     normal distribution will be generated)
%   BVec            - NX1 or 1X1 or 1X2 Ballistic Coefficient of Object (m^2/kg)
%                     (if a 2 element vector is provided, the first element
%                     will be taken as coefficient mean and the second
%                     element will be taken as coefficient variance and a
%                     normal distribution will be generated)
%   SwerlingType    - Text input of Swerling distribution type (optional, Default = 'III')
%                     Allowable Inputs:
%                       * 'I'
%                       * 'II'
%                       * 'III'
%                       * 'IV'
%   Frequency       - 1X1 Frequency of the sensor at which RCS have been
%                     obtained (MHz) (optional, Default = 2000 MHz)
%   NumOfSamples    - 1X1 Number of RCS samples to be generated if only the
%                     median RCS is provided. If Cd and B are also provided
%                     as 1X2 vectors where the first element is the mean
%                     and the second element is the standard deviation,
%                     this will be used to generate normally distributed
%                     vectors for each of them. If the provided RCS is
%                     a Nx1 or the Cd and B vectors are provided as NX1
%                     vectors, number of RCS samples will be set to N
%                     (optional, Default = 1E5)
%
%   Ln_MassCalFac_Mean - 1X1 The mean log-normal mass calibration factor to
%                        scale the calculated masses (optional, Default = 0)
%   Ln_MassCalFac_Sigma- 1X1 The std log-normal mass calibration factor to
%                        scale the calculated masses (optional, Default = 0)
%                        
%   
%
% Outputs:
%   EstimatedMass   - NX1 Estimated Mass sample vector of input object
%   RCS             - NX1 RCS Samples vector used to generate the mass
%                     samples
%   AreaVec         - NX1 Estimated frontal area vector assuming a spherical 
%                     cross section
%
%
% Other m-files required: 	RCSDistribution.m
%                           NASA_SEM_RCSToSizeVec.m
% Subfunctions: None
% MAT-files required: None
%
% See also: none
%
% April 2018; Last revision: 01-Jun-2026
%
% ----------------- BEGIN CODE -----------------
    
    
    % Set Default inputs
    if nargin <3
        error('Insufficient number of inputs')
    elseif nargin == 3
        SwerlingType = [];
        Frequency    = 2000; % Default radar frequency used for estimating object size
        NumOfSamples = 1e5;
        Ln_MassCalFac_Mean = 0;
        Ln_MassCalFac_Sigma = 0;
    elseif nargin == 4
        Frequency    = 2000; % Default radar frequency used for estimating object size
        NumOfSamples = 1e5;
        Ln_MassCalFac_Mean = 0;
        Ln_MassCalFac_Sigma = 0;
    elseif nargin == 5
        NumOfSamples = 1e5;
        Ln_MassCalFac_Mean = 0;
        Ln_MassCalFac_Sigma = 0;
    elseif nargin == 6
        Ln_MassCalFac_Mean = 0;
        Ln_MassCalFac_Sigma = 0;
    end

    if isempty(SwerlingType)
        SwerlingType = 'III';
    end

    % Set Number of samples based on other vectors
    if size(CdVec,1) ~= 1 && size(BVec,1) == size(CdVec,1)
        NumOfSamples = size(CdVec,1);
    end

    % Calculate RCS samples if only one RCS given
    if length(RCS)==1
        RCS = RCSDistribution(RCS,NumOfSamples,SwerlingType);
    end

    % Calculate Cd samples if Standard deviation provided
    if size(CdVec,1)==1 && size(CdVec,2)==2
        Cd = CdVec(1); CdVar = CdVec(2);
        CdVec=abs(Cd+sqrt(CdVar)*randn(NumOfSamples,1)); 
    end

    % Calculate BC samples if Standard deviation provided
    if size(BVec,1)==1 && size(BVec,2)==2
        B = BVec(1); BVar = BVec(2);
        BVec=abs(B+sqrt(BVar)*randn(NumOfSamples,1));
    end
    
    % Set Speed of light
    c = 299792458; % m/s
    
    % radar wavelength in m (associated with RCS frequency); wave equation below
    lambda = c./(Frequency*1e6); 
    
    % Converts RCS value to normalized size value, using the ODPO size estimation model
    x = NASA_SEM_RCSToSizeVec(RCS./lambda.^2);
    
    % Un-normalize the size of the object (meters)
    Size=x.*lambda;

    % Convert Characteristic length to a radius and calculate the frontal
    % area assuming a spherical cross section
    AreaVec=pi*(Size./2).^2;
    
    % vector of satellite mass estimates (from ballistic coefficient equation)
    EstimatedMass = CdVec.*AreaVec./BVec;

    % Scale the estimated masses using the provided calibration factor
    gsamp = Ln_MassCalFac_Mean + Ln_MassCalFac_Sigma*randn(NumOfSamples,1);
    LinMassCalFac = exp(gsamp);
    EstimatedMass = LinMassCalFac .* EstimatedMass;
    

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 04-11-2018 | Initial Development
% S. Es haghi    | 06-01-2026 | Clean out commented segments and add
%                               capability to take sensor frequency and
%                               log-normal mass calibration factors as
%                               input
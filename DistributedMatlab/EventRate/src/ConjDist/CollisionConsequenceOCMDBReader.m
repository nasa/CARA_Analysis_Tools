function [Catastrophic,NumOfPieces,massVec] = CollisionConsequenceOCMDBReader(OCMBDEntry,PrimaryMass,BSigNorm,NumOfSamples)
% CollisionConsequenceOCMDBReader - Calculates which spherical screening 
%                                   volumes are penetrated by a conjuncting 
%                                   secondary
% Syntax: [Catastrophic,NumOfPieces,massVec] = ...
%   CollisionConsequenceOCMDBReader(OCMBDEntry,PrimaryMass,...
%                                   BSigNorm,NumOfSamples)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculates the expected number of resultant debris pieces 
% for an input close approach event for a number of trials enounters. Also 
% outputs the number of these trials with a potentially catastrophic 
% collision as defined in Hejduk et. al. "Consideration of Collision 
% "Consequence" in Satellite Conjunction Assessment and Risk Analysis" 
% 2017.  This is done using an OCMDB Matlab object.
%
% =========================================================================
%
% Input:
%
%   OCMDBEntry      - Single Row of an OCDB object to process
%
%   PrimaryMass     - Mass of Primary Object (kg)
%                     (Optional, default = 2000)
%
%   BSigNorm        - Normalized Ballistic Coefficient Sigma of 
%                     Secondary Object 
%                     (Optional, Default = 0.2 (20%, conservative))
%
%   NumOfSamples    - Number of samples to generate 
%                     (Optional, default = 10000)
%
% =========================================================================
%
% Output:
%
%   Catastrophic    - logical array indicating whether the            [Nx1]
%                     sampled collision is catastrophic
%
%   NumOfPieces     - array of the number of pieces expected to       [Nx1]
%                     be generated from a collision for each sample
%
%   massVec         - array of the secondary object mass              [Nx1]
%                     estimates for each individual sample
%
% =========================================================================
%
% References:
%
%   Hejduk et. al. "Consideration of Collision "Consequence" in Satellite 
%   Conjunction Assessment and Risk Analysis" 2017
%
% =========================================================================
%
% Dependencies:
%
%   RCSDistribution.m
%   EstimateMassFromRCS.m
%   CollisionConsequenceNumPieces.m
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Feb 2018;  Latest update: Nov 2020
%
% ----------------- BEGIN CODE -----------------
    persistent pathsAdded
    if isempty(pathsAdded)
        % Paths for collision consequence mode
        s = what(fullfile(p,'../Utils/AugmentedMath')); addpath(s.path);
        s = what(fullfile(p,'../ProbabilityOfCollision')); addpath(s.path);
        pathsAdded = true;
    end
    %% Set up Defaults
    % Default number of samples
    if nargin < 4
        NumOfSamples = [];
    end
    % Default sigma(B)/B value
    if nargin < 3 || isempty(BSigNorm)
        BSigNorm = 0.2; 
    end
    % Default primary object mass at collision time
    if nargin < 2 || isempty(PrimaryMass)
        PrimaryMass = 2000;
    end
    
    % Set Cd Values
    DefaultCd       = 2.7; % Default Drag Coefficient Used for Estimating Masses of Primary and Secondary Objects
    DefaultCdSigma  = 0.10; % Relative Cd uncertainty for use in Mass Estimation Process
    Cd          = DefaultCd;
    CdVar       = (DefaultCd*DefaultCdSigma)^2; % Get Cd Variance from Relative Uncertainty Specified above
    
    % Set B value for input
    B = OCMBDEntry(1,116);
    
    % Set negative B values to zero
    B = max(0,B);
    
    % Set BVar to be the estimated variance of B = BSig^2
    BVar = (BSigNorm*B)^2;
    
    % Get relative velocity magnitude
    % VRel = norm(OCMBDEntry(1,24:26));
    VRel = OCMBDEntry(1,20);    
    
    % Set RCS Value
    RCS = OCMBDEntry(1,114);

    % Estimate the 99% confidence upper-limit secondary mass 
    [massVec,~] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,0.99,NumOfSamples);
    
    % Estimate the catastrophic nature and the number of pieces for the
    % primary-secondary collision, using the upper-limit secondary mass
    [Catastrophic,NumOfPieces] = CollisionConsequenceNumPieces(PrimaryMass,VRel,massVec);
    
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 02-05-2018 | Initial Development
% D. Hall        | 02-14-2019 | Set negative OCMBD B values to zero, and
%                               made BVar hold variance by squaring BSig
% D. Hall        | 11-04-2020 | Worked with T. Lechtenberg to interface
%                               properly with latest versions of functions
%                               EstimateObjectMassQuantiles.m and
%                               CollisionConsequenceNumPieces.m. Changed
%                               the default NumOfSamples to be set within 
%                               function EstimateObjectMassQuantiles.m.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
function [Catastrophic,NumOfPieces,MeanNumOfPieces,out] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass,Lc)
%
% CollisionConsequenceNumPieces - Calculates the expected number of
% resultant debris pieces for an input close approach event for a number of
% trials enounters.  Also outputs the number of these trials with a
% potentially catastrophic collision as defined in Hejduk et. al.
% "Consideration of Collision "Consequence" in Satellite Conjunction
% Assessment and Risk Analysis" 2017
%
% Syntax:   [Catastrophic,NumOfPieces,ExpectedNumOfPieces,out] = ...
%           CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass,Lc)
%
% =========================================================================
%
% Copyright (c) 2019-2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Inputs:
%
%   PrimaryMass     - 1X1 or NX1 Mass of Primary Object (kg) (The primary
%                     mass input can be a vector the same size as the
%                     secondary mass vector)
%   VRel            - 1X1, 3X1, 1X3 or NX1 vector of the relative velocity
%                     between the primary and secondary objects (m/s)
%                     (relative velocity uncertainty can be taken into
%                     account if an N by 1 relative velocity vector is
%                     provided as input)
%   SecondaryMass   - 1X1 or NX1 Matrix of Secondary Mass Values for Examination
%                     (kg)
%   Lc              - Characteristic Length of Debris Pieces to use as
%                     threshold for reporting (m) (Optional,defaults to
%                     0.05 (5 cm))
%
% =========================================================================
%
% Outputs:
%
%   Catastrophic    - [NX1] logical array indicating whether the
%                     sampled collision is catastrophic
%   NumOfPieces     - [NX1] array of the number of pieces
%                     to be generated from a collision for each sample
%                     (rounded numbers)
%   MeanNumOfPieces - [1X1] statistically expected number of fragments,
%                     given a direct-hit collision
%   out             - output structure with auxiliary information:
%                      * NumOfPiecesGivenCollStDev: NumOfPieces standard
%                        deviation
%                      * NumOfPiecesGivenCollStDevMean: NumOfPieces standard
%                        deviation of the mean given the number of samples
%
% =========================================================================
%
% Other m-files required: 	None
% Subfunctions:             None
% MAT-files required:       None
%
% See also: Hall, D. and Baars, L., "Satellite Collision and Fragmentation
%           Probabilities Using Radar-Based Size and Mass Estimates,"
%           Journal of Spacecraft and Rockets, 2023
% 
%           Lechtenberg, T., "An Operational Algorithm for Evaluating
%           Satellite Collision Consequence," AAS Astrodynamics Specialist
%           Conference, 2019, AAS 19-669
%
%           Krisko, P., "Proper Implementation of the 1998 NASA Breakup
%           Model," Orbital Debris Quarterly News, Volume 15, Issue 4,
%           October 2011.
%
% =========================================================================
%
% Initial version: May 2019;  Latest update: June 2026
%
% ----------------- BEGIN CODE -----------------
    
    % Set up Defaults
    if nargin < 4 || isempty(Lc) || Lc==0
        Lc = 0.05;
    end
    
    out = [];
    
    % Set Default inputs
    CatastrophicThreshold   = 40000; % Threshold for Catastrophic Collision (Joules/kg)
    
    % Get relative velocity magnitude
    if min(size(VRel))==1 && max(size(VRel))==3
        VRel                = norm(VRel);
    else
        VRel                = abs(VRel);
    end

    % Number of Samples
    NoS = max([numel(VRel),numel(SecondaryMass),numel(PrimaryMass)]);

    % Reshape Inputs
    VRel                    = reshape(VRel,numel(VRel),1);
    PrimaryMass             = reshape(PrimaryMass,numel(PrimaryMass),1);
    SecondaryMass           = reshape(SecondaryMass,numel(SecondaryMass),1);

    % Make primary mass vector if a single element provided
    if isscalar(PrimaryMass) && NoS~=numel(PrimaryMass)
        PrimaryMass = repmat(PrimaryMass,[NoS 1]);
    end

    % Make secondary mass vector if a single element provided
    if isscalar(SecondaryMass) && NoS~=numel(SecondaryMass)
        SecondaryMass = repmat(SecondaryMass,[NoS 1]);
    end

    % Make relative velocity vector if only one relative velocity is given
    if isscalar(VRel) && numel(VRel)~=NoS
        VRel = repmat(VRel,[NoS 1]);
    end

    % Check if the dimension of all inputs agree with eachother
    if numel(VRel)~=numel(SecondaryMass) || numel(VRel)~=numel(PrimaryMass)
        error('Input dimensions mismatch');
    end

    % ODPO relative velocity kinetic energy equation. (adjusted to use larger mass as dividend) 
    CollisionEnergy            = zeros(size(SecondaryMass)); %Preallocate
    ind_PbS                    = SecondaryMass<=PrimaryMass;
    CollisionEnergy(ind_PbS)   = 0.5.*SecondaryMass(ind_PbS).*VRel(ind_PbS).^2./PrimaryMass(ind_PbS);
    CollisionEnergy(~ind_PbS)  = 0.5.*PrimaryMass(~ind_PbS).*VRel(~ind_PbS).^2./SecondaryMass(~ind_PbS);
    
    % catastrophic / non-catastrophic determination
    Catastrophic            = CollisionEnergy>CatastrophicThreshold;
    
    % in a catastrophic collision, BigM is defined as the sum of the two objects' masses; in a
    % non-catastrophic collision, BigM is defined as the product of the mass of the smaller object and
    % the collision velocity (in km/sec)
    BigM=Catastrophic.*(SecondaryMass+PrimaryMass)+~Catastrophic.*min(SecondaryMass,PrimaryMass).*(VRel/1000).^2;
    % formula for number of pieces from ODPO model
    NumOfPieces=0.1.*BigM.^0.75.*(Lc).^-1.71;

    % Calculate statistically expected number of fragments in case of
    % collision
    MeanNumOfPieces = mean(NumOfPieces);
    NumOfPiecesGivenCollStDev = std(NumOfPieces);
    NumOfPiecesGivenCollStDevMean = NumOfPiecesGivenCollStDev/sqrt(numel(NumOfPieces));
    out.Lc = Lc;
    out.CatastrophicEnergyThreshold = CatastrophicThreshold;
    out.NumOfPiecesGivenCollStDev = NumOfPiecesGivenCollStDev;
    out.NumOfPiecesGivenCollStDevMean = NumOfPiecesGivenCollStDevMean;

    % Round the vector of number of fragments
    NumOfPieces = round(NumOfPieces);

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 05-23-2019 | Initial Development
% T. Lechtenberg | 10-19-2021 | Addition of reported Debris Piece sizes as
%                               independent variable
% L. Baars       | 09-26-2025 | Fixed calculation of BigM per ODQN 15-4
%                               (corrections to NASA breakup model,
%                               equation 4)
% S. Es haghi    | 05-29-2026 | Output statistically expected number of pieces
%                               and ability to accept Primary mass vector
% S. Es haghi    | 06-08-2026 | Account for velocity vector input (velocity
%                               uncertainty included) and update header
%
% =========================================================================
%
% Copyright (c) 2019-2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

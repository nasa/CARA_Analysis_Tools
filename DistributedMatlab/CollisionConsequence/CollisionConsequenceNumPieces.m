function [Catastrophic,NumOfPieces] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass,Lc)
%
% CollisionConsequenceNumPieces - Calculates the expected number of
% resultant debris pieces for an input close approach event for a number of
% trials enounters.  Also outputs the number of these trials with a
% potentially catastrophic collision as defined in Hejduk et. al.
% "Consideration of Collision "Consequence" in Satellite Conjunction
% Assessment and Risk Analysis" 2017
%
% Syntax:   [Catastrophic,NumOfPieces] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass,Lc)
%
% =========================================================================
%
% Copyright (c) 2019-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Inputs:
%
%   PrimaryMass     - 1X1 Mass of Primary Object (kg)
%   VRel            - 1x1, 3X1, or 1X3 vector of the relative velocity
%                     between the primary and secondary objects (m/s)
%   SecondaryMass   - NX1 Matrix of Secondary Mass Values for Examination
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
%                     expected to be generated from a collision for each sample
%
% =========================================================================
%
% Other m-files required: 	None
% Subfunctions:             None
% MAT-files required:       None
%
% See also: Lechtenberg, T., "An Operational Algorithm for Evaluating
%           Satellite Collision Consequence," AAS Astrodynamics Specialist
%           Conference, 2019, AAS 19-669
%
%           Krisko, P., "Proper Implementation of the 1998 NASA Breakup
%           Model," Orbital Debris Quarterly News, Volume 15, Issue 4,
%           October 2011.
%
% =========================================================================
%
% Initial version: May 2019;  Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------
    
    %% Set up Defaults
    if nargin < 4 || isempty(Lc) || Lc==0
        Lc = 0.05;
    end
    
    % Set Default inputs
    CatastrophicThreshold   = 40000; % Threshold for Catastrophic Collision (Joules/kg)
    
    % Get relative velocity magnitude
    VRel                    = norm(VRel);
    
    % Reshape Secondary Mass Input
    SecondaryMass           = reshape(SecondaryMass,numel(SecondaryMass),1);
    
    % ODPO relative velocity kinetic energy equation. (adjusted to use larger mass as dividend) 
    CollisionEnergy                             = zeros(size(SecondaryMass)); %Preallocate
    CollisionEnergy(SecondaryMass<=PrimaryMass) = 0.5.*SecondaryMass(SecondaryMass<=PrimaryMass).*VRel.^2./PrimaryMass;
    CollisionEnergy(SecondaryMass>PrimaryMass)  = 0.5.*PrimaryMass.*VRel.^2./SecondaryMass(SecondaryMass>PrimaryMass);
    
    % catastrophic / non-catastrophic determination
    Catastrophic            = CollisionEnergy>CatastrophicThreshold;
    
    % in a catastrophic collision, BigM is defined as the sum of the two objects' masses; in a
    % non-catastrophic collision, BigM is defined as the product of the mass of the smaller object and
    % the collision velocity (in km/sec)
    BigM=Catastrophic.*(SecondaryMass+PrimaryMass)+~Catastrophic.*min(SecondaryMass,PrimaryMass).*(VRel/1000)^2;
    % formula for number of pieces from ODPO model
    NumOfPieces=round(0.1.*BigM.^0.75.*(Lc).^-1.71);

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

% =========================================================================
%
% Copyright (c) 2019-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

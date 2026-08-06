function J = JacobianEphoffset2CartesianFull(t,Et,MeanXt,eph)
% JacobianEphoffset2CartesianFull - Calculate the Jacobian, J = dX/dE,
%                                   between the ephemeris-offset elements
%                                   and the cartesian state.
%
% Syntax: J = JacobianEphoffset2CartesianFull(t, Et, MeanXt, eph);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    t          -   Ephemeris time
%
%    Et         -   Ephemeris-offset (EO) state at time t   [6x1]
%                   (km and km/s)
%
%    MeanXt     -   Mean cartesian state at time t          [6x1]
%                   (km and km/s)
%
%    eph        -   Ephemeris table structure with fields:
%                     .T - Ephemeris times                  [1xN]
%                     .X - State vectors                    [6xN]
%                     .N - Number of points
%
% =========================================================================
%
% Output:
%
%   J           -   Jacobian matrix dX/dE                   [6x6]
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Initialize Jacobian Matrix
J = NaN(6,6);

% Perturbations km and km/s
epsE = [1e-4,1e-4,1e-4,1e-7,1e-7,1e-7];
twoepsE = 2*epsE;

% Numerically compute Jacobian Matrix using central differences method
for i = 1:6

    % Re-set/initialize perturbed states
    Eminus = Et;
    Eplus  = Et;

    % Perturb state
    Eminus(i) = Et(i) - epsE(i);
    Eplus (i) = Et(i) + epsE(i);

    % Convert the perturbed Cartesian states into Equinoctial states 
    Xminus = convert_ephoffset_to_cartesian(t,Eminus,MeanXt,eph);
    Xplus  = convert_ephoffset_to_cartesian(t,Eplus ,MeanXt,eph);

    % Partial derivatives based on central differences
    J(:,i) = (Xplus(1:6) - Xminus(1:6)) / twoepsE(i);

end
    
return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
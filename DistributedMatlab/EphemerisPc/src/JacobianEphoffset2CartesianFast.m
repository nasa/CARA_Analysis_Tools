function J = JacobianEphoffset2CartesianFast(t,Et,AuxEt,MeanXt,eph)
% JacobianEphoffset2CartesianFast - Calculate the Jacobian, J = dX/dE,
%                                   between the ephemeris-offset elements
%                                   and the cartesian state.
%
% Syntax: J = JacobianEphoffset2CartesianFast(t, Et, AuxEt, MeanXt, eph);
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
%    AuxEt      -   Auxiliary Et quantities (i.e., the 'out'
%                   structure produced by the function
%                   convert_cartesian_to_ephoffset.m)
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
J = zeros(6,6);

% Velocity-velocity part of Jacobian matrix 
J(4:6,4:6) = [AuxEt.Xhat_tau; ...
              AuxEt.Yhat_tau; ...
              AuxEt.Zhat_tau]';

% Perturbations km and km/s
epsE13 = 1e-4; twoepsE13 = 2*epsE13;

% Numerically compute Jacobian Matrix using central differences method
for i = 1:3

    % Re-set/initialize perturbed states
    Eminus = Et;
    Eplus  = Et;

    % Perturb state
    Eminus(i) = Et(i) - epsE13;
    Eplus (i) = Et(i) + epsE13;

    % Convert the perturbed Cartesian states into Equinoctial states 
    Xminus = convert_ephoffset_to_cartesian(t,Eminus,MeanXt,eph);
    Xplus  = convert_ephoffset_to_cartesian(t,Eplus ,MeanXt,eph);

    % Partial derivatives based on central differences
    J(:,i) = (Xplus(1:6) - Xminus(1:6)) / twoepsE13;

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
function J = JacobianCartesian2EphoffsetFast(t,Xt,AuxEt,MeanXt,eph)
%
% Calculate the Jacobian, J = dE/dX, between the ephemeris-offset elements
% and the cartesian state
%
% INPUT:
% 
%   t  = Ephemeris time
%   Xt = Cartesian state at time t; Xt = (rt',vt')' [ km & km/s ]
%   Et = Ephemeris-offset elements at time t        [ km & km/s ]
%   AuxEt = Auxiliary Et quantities
%           (i.e., the 'out' structure produced by the function
%           convert_cartesian_to_ephoffset.m)
%   MeanXt = mean Cartesian state at time t         [ km & km/s ]
%   eph = Ephemeris table structure
%
% OUTPUT:
%
%   J = dE/dX [6x6] = Jacobian matrix
%

% Persistent variable for warning
persistent warning_issued
if isempty(warning_issued)
    warning_string = ['Calculating JacobianCartesian2EphoffsetFast is ' ...
        'usually less efficient than calculating the inverse of ' ...
        'JacobianEphoffset2CartesianFast'];
    warning(warning_string);
    warning_issued = true;
end

% Initialize Jacobian Matrix
J = zeros(6,6);

% Ensure that AuxEt is not empty
if isempty(AuxEt)
    [~,AuxEt] = convert_cartesian_to_ephoffset(t,Xt,MeanXt,eph);
end

% Velocity-velocity part of Jacobian matrix 
J(4:6,4:6) = [AuxEt.Xhat_tau; ...
              AuxEt.Yhat_tau; ...
              AuxEt.Zhat_tau]';

% Perturbations for first three eph.offset elements, E(1:3), in km
epsX13 = 1e-4; twoepsX13 = 2*epsX13;

% Numerically calculate the required derivatives
for i = 1:3

    % Re-set/initialize perturbed states
    Xminus = Xt;
    Xplus  = Xt;

    % Perturb state
    Xminus(i) = Xt(i) - epsX13;
    Xplus (i) = Xt(i) + epsX13;

    % Convert the perturbed Cartesian states into Equinoctial states 
    [Eminus,~] = convert_cartesian_to_ephoffset(t,Xminus,MeanXt,eph);
    [Eplus ,~] = convert_cartesian_to_ephoffset(t,Xplus ,MeanXt,eph);

    % Partial derivatives based on central differences
    J(:,i) = (Eplus(1:6) - Eminus(1:6)) / twoepsX13;

end
    
return
end
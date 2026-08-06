function J = JacobianCartesian2EphoffsetFull(t,Xt,MeanXt,eph)
%
% Calculate the Jacobian, J = dE/dX, between the ephemeris-offset elements
% and the cartesian state
%
% Note: It is likely more efficient to calculate J = (dX/dE)^-1, using the
% function JacobianEphoffset2CartesianFull to calculate dX/dE. The reason
% for this is that the latter function needs to interpolate from the
% ephemeris table much less frequently
%
% INPUT:
% 
%   t  = Ephemeris time
%   Xt = Cartesian state at time t = (r',v')' [ km & km/s ]
%   MeanXt = mean Cartesian state at time t   [ km & km/s ]
%   eph = Ephemeris table structure
%
% OUTPUT:
%
%   J = dE/dX [6x6] = Jacobian matrix
%

% Persistent variable for warning
persistent warning_issued
if isempty(warning_issued)
    warning_string = ['Using JacobianCartesian2EphoffsetFull is ' ...
        'usually less efficient than using the inverse of ' ...
        'JacobianEphoffset2CartesianFull'];
    warning(warning_string);
    warning_issued = true;
end

% Initialize Jacobian Matrix
J = NaN(6,6);

% Perturbations km and km/s
epsX = [1e-4,1e-4,1e-4,1e-7,1e-7,1e-7];
twoepsX = 2*epsX;

% Numerically compute Jacobian Matrix using central differences method
for i = 1:6

    % Re-set/initialize perturbed states
    Xminus = Xt;
    Xplus  = Xt;

    % Perturb state
    Xminus(i) = Xt(i) - epsX(i);
    Xplus (i) = Xt(i) + epsX(i);

    % Convert the perturbed Cartesian states into Equinoctial states 
    Eminus = convert_cartesian_to_ephoffset(t,Xminus,MeanXt,eph);
    Eplus  = convert_cartesian_to_ephoffset(t,Xplus ,MeanXt,eph);

    % Partial derivatives based on central differences
    J(:,i) = (Eplus(1:6) - Eminus(1:6)) / twoepsX(i);

end
    
return
end
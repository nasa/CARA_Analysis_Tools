function [EQ] = Cart2Eq(CART)
%
% Cart2Eq  - Converts an object's Cartesian elements to Equinoctial elements.
%            The state vector must be input in units of km and km/s. These
%            equinoctial elements follow the same convention used on the VCM 
%            equinoctial covariance matrix.
%            The conversion is valid for circular, elliptic, and hyperbolic
%            orbits. These equinoctial elements exhibit no singularity for 
%            zero eccentricity and orbital inclinations equal to 0 and 90 
%            degrees. However, two of the components (h and k approach 
%            +/- infinity) are singular for an orbit inclination of 180 degrees.
%
% Synthax:  [f,g,L,n,k,h] = Cart2Eq(r,v)
%
% Inputs: 
%    r     - Position vector in ECI coordinates.          Units: km 
%            (1x3 row vector or 3x1 column vector)
%    v     - Velocity vector in ECI coordinates.          Units: km/s
%            (1x3 row vector or 3x1 column vector)
%
% Outputs:
%    f     - e*cos(w+Omega)                               Units: None
%    g     - e*sin(w+Omega)                               Units: None
%    L     - Omega+w+M (mean longitude)                   Units: rad   [0, 2pi]
%    n     - sqrt(mu/a^3) (mean motion)                   Units: 1/s
%    k     - tan(i/2)*sin(Omega)                          Units: None
%    h     - tan(i/2)*cos(Omega)                          Units: None
%
% Examples/Validation Cases:
%
% Other m-files required: Cart2Kep.m, Kep2Eq.m
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% March 2013; Last Revision: 11-Dec-2014
%
% ----------------- BEGIN CODE -----------------
    
    % Convert Cartesian elements to Keplerian elements
    [KEP] = Cart2Kep(CART,'Mean','Rad');
    %[a,e,i,Omega,w,M] = Cart2Kep(r,v,'Mean','Rad');
    
    % Convert Keplerian elements to Equinoctial elements
    [EQ]     = Kep2Eq(KEP);

return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% D. Plakalovic  | Mar - 2013 |  Initial Development
% D. Plakalovic  | 12-11-2014 |  Added more detailed commentary
%
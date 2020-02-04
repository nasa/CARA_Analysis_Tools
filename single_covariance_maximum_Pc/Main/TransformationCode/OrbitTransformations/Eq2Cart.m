function [CART] = Eq2Cart(EQ)

% EQ   is nx6 matrix of Equinoctial states             (by row)
% KEP  is nx6 matrix of Keplerian corresponding states (by row)
% CART is nx6 matrix of Cartesian corresponding state  (by row)


% function [r,v] = Eq2Cart(f,g,L,n,k,h)
%
% Eq2Cart - Converts an object's Equinoctial elements to Cartesian elements.
%           These equinoctial elements follow the same convention used on 
%           the VCM equinoctial covariance matrix. Note that L is mean
%           longitude with respect to mean anomaly.
%           The transformation is valid for circular, elliptic, and
%           hyperbolic orbits. These equinoctial elements exhibit no 
%           singularity for zero eccentricity and orbital inclinations equal
%           to 0 and 90 degrees. However, two of the components are singular
%           for an orbit inclination of 180 degrees (h and k approach +/- infinity).
%
% Synthax:  [r,v] = Eq2Cart(f,g,L,n,k,h)
%
% Inputs: 
%    f     - e*cos(w+Omega)                                Units: None
%    g     - e*sin(w+Omega)                                Units: None
%    L     - Omega+w+M (mean longitude)                    Units: rad   [-2pi,2pi]
%    n     - sqrt(mu/a^3) (mean motion)                    Units: 1/s
%    k     - tan(i/2)*sin(Omega)                           Units: None
%    h     - tan(i/2)*cos(Omega)                           Units: None
%
% Outputs:
%    r - Position vector in ECI coordinates           Units: km 
%        (1x3 row vector)
%    v - Velocity vector in ECI coordinates           Units: km/s
%        (1x3 row vector)
%
% Examples/Validation Cases:
%
% Other m-files required: Eq2Kep.m, Kep2Cart.m
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% March 2013; Last Revision: 11-Dec-2014
%
% ----------------- BEGIN CODE -----------------

    % Convert Equinoctial elements to Keplerian elements
    [KEP] = Eq2Kep(EQ);
    % [a,e,i,Omega,w,M] = Eq2Kep(f,g,L,n,k,h);

    % Convert Keplerian elements to Cartesian elements
    [CART] = Kep2Cart(KEP);
    % [r,v] = Kep2Cart(a,e,i,Omega,w,M);

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
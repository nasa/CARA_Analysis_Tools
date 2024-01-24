function [KEP] = Eq2Kep(EQ)
% Eq2Kep - Converts an object's Equinoctial elements to Keplerian elements
%
% Syntax: [a,e,i,Omega,w,M] = Eq2Kep(f,g,L,n,k,h);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Converts an object's Equinoctial elements to Keplerian elements. These 
%   equinoctial elements follow the same convention used on the VCM 
%   equinoctial covariance matrix. Note that L is mean longitude with 
%   respect to mean anomaly. The transformation is valid for circular, 
%   elliptic, and hyperbolic orbits. These equinoctial elements exhibit no 
%   singularity for zero eccentricity and orbital inclinations equal to 0 
%   and 90 degrees. However, two of the components are singular for an 
%   orbit inclination of 180 degrees (h and k approach +/- infinity).
%
% =========================================================================
%
% Input:
%
%    f      -   e*cos(w+Omega)                                
%    g      -   e*sin(w+Omega)                               
%    L      -   Omega+w+M (mean longitude) (rad)
%    n      -   sqrt(mu/a^3) (mean motion) (1/s)
%    k      -   tan(i/2)*sin(Omega)
%    h      -   tan(i/2)*cos(Omega)
%
% =========================================================================
%
% Output:
%
%    a      -   Semi-Major axis (km)
%    e      -   Eccentricity                                  
%    i      -   Inclination (rad)
%    Omega  -   Right Ascension of the Ascending Node (RAAN) (rad)
%    w      -   Argument of Perigee (rad)
%    M      -   Mean Anomaly (rad) 
%
% =========================================================================
%
% Initial version: Mar 2013;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

    % EQ = [f,g,L,n,k,h,fr]
    
    f = EQ(:,1);
    g = EQ(:,2);
    L = EQ(:,3);
    n = EQ(:,4);
    k = EQ(:,5);
    h = EQ(:,6);
    
    if size(EQ,2) < 7
        fr = 1;
    else
        fr = EQ(:,7);
    end

    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu = 3.986004418e5;

    % Semi-Major axis
    a = (mu./n.^2).^(1/3);
    
    % Eccentricity
    e = sqrt(f.^2+g.^2);
    
    if any(abs(e - 1) <= 1E-5)
        error('Eq2Kep:ParabolicOrRectilinearOrbit', 'Parabolic or rectilinear orbit - cannot calculate equinoctial elements');
    end
    
    % Inclination
%     i = atan2(2*sqrt(h.^2 + k.^2),1-h.^2-k.^2);
    i = pi*(1-fr)/2+2.*fr.*atan(sqrt(h.^2 + k.^2));
    
    % Right ascension of the ascending node (RAAN)
    Omega = atan2(k,h);
    
    Omega = mod(Omega, 2 * pi);
    id = Omega >= -pi & Omega < 0;
    Omega(id) = Omega(id) + 2*pi;
        
    % Argument of perigee
%     w = atan2(g.*h-f.*k,f.*h+g.*k);
    w = atan2(g,f) - fr.*atan2(k,h);
    
    w = mod(w, 2 * pi);
    id = w >= -pi & w < 0;
    w(id) = w(id) + 2*pi;
        
    % Mean anomaly
    M = mod(L - fr.*Omega - w, 2*pi);
    
    KEP = [a,e,i,Omega,w,M];

return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% D. Plakalovic  | Mar - 2013 | Initial Development
% D. Plakalovic  | 04-02-2013 | Added a mod function to L-component in
%                               order to properly handle angles greater
%                               than 360 degrees. Also added if-statements
%                               to 'Omega' and 'w' components to account
%                               for negative angles, so that they display
%                               0 to 360 degrees
% D. Plakalovic  | 12-11-2014 | Added more detailed commentary
% E. White       | 06-30-2023 | Updated documentation, removed unused
%                               commented-out code, added check for
%                               parabolic and rectilinear orbits

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
    
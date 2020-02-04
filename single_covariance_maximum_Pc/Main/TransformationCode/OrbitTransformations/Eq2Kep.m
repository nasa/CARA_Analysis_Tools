function [KEP] = Eq2Kep(EQ)
% function [a,e,i,Omega,w,M] = Eq2Kep(f,g,L,n,k,h)
%
% Eq2Kep   - Converts an object's Equinoctial elements to Keplerian elements.
%            These equinoctial elements follow the same convention used on 
%            the VCM equinoctial covariance matrix. Note that L is mean
%            longitude with respect to mean anomaly.
%            The transformation is valid for circular, elliptic, and
%            hyperbolic orbits. These equinoctial elements exhibit no 
%            singularity for zero eccentricity and orbital inclinations equal
%            to 0 and 90 degrees. However, two of the components are singular
%            for an orbit inclination of 180 degrees (h and k approach +/- infinity).
%
% Synthax:   [a,e,i,Omega,w,M] = Eq2Kep(f,g,L,n,k,h)
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
%    a     - Semi-Major axis                               Units: km
%    e     - Eccentricity                                  Units: None
%    i     - Inclination                                   Units: rad   [0,pi]
%    Omega - Right Ascension of the Ascending Node (RAAN)  Units: rad   [0,2pi] 
%    w     - Argument of Perigee                           Units: rad   [0,2pi]
%    M     - Mean Anomaly                                  Units: rad   [0,2pi]
%
% Examples/Validation Cases:
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% March 2013; Last Revision: 11-Dec-2014
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
    
    % Inclination
%     i = atan2(2*sqrt(h.^2 + k.^2),1-h.^2-k.^2);
    i = pi*(1-fr)/2+2.*fr.*atan(sqrt(h.^2 + k.^2));
    
    % Right ascension of the ascending node (RAAN)
    Omega = atan2(k,h);
    
    id = Omega >= -pi & Omega < 0;
    Omega(id) = Omega(id) + 2*pi;
    
%     if (Omega >= -pi && Omega < 0)
%         Omega = Omega + 2*pi;
%     end
        
    % Argument of perigee
%     w = atan2(g.*h-f.*k,f.*h+g.*k);
    w = atan2(g,f) - fr.*atan2(k,h);
    
    id = w >= -pi & w < 0;
    w(id) = w(id) + 2*pi;
    
%     if (w >= -pi && w < 0)
%         w = w + 2*pi;
%     end
        
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
% D. Plakalovic  | Mar - 2013 |  Initial Development
% D. Plakalovic  | 04-02-2013 |  Added a mod function to L-component in
%                                order to properly handle angles greater
%                                than 360 degrees. Also added if-statements
%                                to 'Omega' and 'w' components to account
%                                for negative angles, so that they display
%                                0 to 360 degrees
% D. Plakalovic  | 12-11-2014 |  Added more detailed commentary
%
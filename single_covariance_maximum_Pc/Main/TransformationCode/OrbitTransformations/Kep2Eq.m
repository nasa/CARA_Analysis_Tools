function [EQ] = Kep2Eq(KEP)
%
% Kep2Eq   - Converts an object's Keplerian elements to Equinoctial elements.
%            These equinoctial elements follow the same convention used on 
%            the VCM equinoctial covariance matrix. Note that the input
%            anomaly angle must be with respect to mean anomaly, and all 
%            angles inputs are in the units of radians.
%            The transformation is valid for circular, elliptic, and
%            hyperbolic orbits. These equinoctial elements exhibit no 
%            singularity for zero eccentricity and orbital inclinations 
%            equal to 0 and 90 degrees. However, two of the components are 
%            singular for an orbit inclination of 180 degrees 
%            (h and k approach +/- infinity).
%
% Synthax:   [f,g,L,n,k,h] = Kep2Eq(a,e,i,Omega,w,M)
%
% Inputs: 
%    a     - Semi-Major axis                              Units: km
%    e     - Eccentricity                                 Units: None
%    i     - Inclination                                  Units: rad   [-2pi, 2pi]
%    Omega - Right Ascension of the Ascending Node (RAAN) Units: rad   [-2pi, 2pi]
%    w     - Argument of Perigee                          Units: rad   [-2pi, 2pi]
%    M     - Mean Anomaly                                 Units: rad   [-2pi, 2pi]
%
% Outputs:
%    f     - af  = e*cos(w+Omega)                               Units: None
%    g     - ag  = e*sin(w+Omega)                               Units: None
%    L     - Omega+w+nu (true longitude)                        Units: rad   [0, 2pi]
%    n     - sqrt(mu/a^3) (mean motion)                         Units: 1/s
%    k     - chi = tan(i/2)*sin(Omega)                          Units: None
%    h     - psi = tan(i/2)*cos(Omega)                          Units: None
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

    a     = KEP(:,1);
    e     = KEP(:,2);
    i     = KEP(:,3);
    Omega = KEP(:,4);
    w     = KEP(:,5);
    M     = KEP(:,6);

    % Get Retrograde Factor
    fr = 1;
    if abs(i) > 0.95*pi && abs(i)<1.05*pi
        fr = -1;
    end
    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu = 3.986004418e5;

    % f - component
    f = e * cos(w+fr*Omega);
    
    % g - component
    g = e * sin(w+fr*Omega);
    
    % True Longitude
    L = mod(fr*Omega + w + M, 2*pi);
    
    % Mean motion
    n  = sqrt(mu/a^3);
    
    % k - component
    if fr == -1
        k = atan(i/2) * sin(Omega);
    else
        k = tan(i/2) * sin(Omega);
    end
    k = tan(i/2).^fr * sin(Omega);
    
    % h - component
    if fr == -1
        h = atan(i/2) * cos(Omega);
    else
        h = tan(i/2) * cos(Omega);
    end
    h = tan(i/2).^fr * cos(Omega);
    
    EQ = [f,g,L,n,k,h,fr];
    
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
%                                than 360 degrees
% D.Plakalovic   | 12-11-2014 |  Modified function to streamline with
%                                Astrostd VCM covariance convention
%
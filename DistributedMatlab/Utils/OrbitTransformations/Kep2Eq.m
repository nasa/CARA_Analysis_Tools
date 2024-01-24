function [EQ] = Kep2Eq(KEP)
% Kep2Eq - Converts an object's Keplerian elements to Equinoctial elements
%
% Syntax: [f,g,L,n,k,h,fr] = Kep2Eq(a,e,i,Omega,w,M);
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
%   Converts an object's Keplerian elements to Equinoctial elements. These 
%   equinoctial elements follow the same convention used on the VCM 
%   equinoctial covariance matrix. Note that the input anomaly angle must 
%   be with respect to mean anomaly, and all angles inputs are in the units
%   of radians. The transformation is valid for circular, elliptic, and
%   hyperbolic orbits. These equinoctial elements exhibit no singularity 
%   for zero eccentricity and orbital inclinations  equal to 0 and 90 
%   degrees. However, two of the components are singular for an orbit 
%   inclination of 180 degrees  (h and k approach +/- infinity).
%
% =========================================================================
%
% Input:
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
% Output:
%
%    f      -   af  = e*cos(w+Omega)
%    g      -   ag  = e*sin(w+Omega)
%    L      -   Omega+w+M (mean longitude) (rad)          
%    n      -   sqrt(mu/a^3) (mean motion) (1/s)
%    k      -   chi = tan(i/2)*sin(Omega)  
%    h      -   psi = tan(i/2)*cos(Omega)
%    f      -   retrograde factor (1 or -1)
%
% =========================================================================
%
% Initial version: Mar 2013;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

    a     = KEP(1);
    e     = KEP(2);
    i     = KEP(3);
    Omega = KEP(4);
    w     = KEP(5);
    M     = KEP(6);
    
    if abs(e - 1) <= 1E-5
        error('Kep2Eq:ParabolicOrRectilinearOrbit', 'Parabolic or rectilinear orbit - cannot calculate equinoctial elements');
    end
    
    if e > 1
        error('Kep2Eq:HyperbolicOrbit', 'Hyperbolic orbit - cannot calculate equinoctial elements');
    end

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
    k = tan(i/2).^fr * sin(Omega);
    
    % h - component
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
% D. Plakalovic  | Mar - 2013 | Initial Development
% D. Plakalovic  | 04-02-2013 | Added a mod function to L-component in
%                               order to properly handle angles greater
%                               than 360 degrees
% D.Plakalovic   | 12-11-2014 | Modified function to streamline with
%                               Astrostd VCM covariance convention
% E. White       | 06-30-2023 | Added compliant documentation, removed
%                               broken vectorization, removed unused and
%                               incorrect code blocks

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

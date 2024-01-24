function [CART] = Kep2Cart(KEP)
% Kep2Cart - Converts an object's Keplerian elements to Cartesian elements
%
% Syntax: [CART] = Kep2Cart(a,e,i,Omega,w,M);
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
%   Converts an object's Keplerian elements to Cartesian elements. Note 
%   that the input anomaly angle must be with respect to mean anomaly and 
%   all angle inputs are in the units of radians.
%
%   NOTE: Mean2EccAnomaly.m is a dependency for this function.
%
% =========================================================================
%
% Input:
%
%   a       -   Semi-major axis (km)                                  [nx1]
%   e       -   Eccentricity                                          [nx1]
%   i       -   Inclination (rad)                                     [nx1]
%   Omega   -   Right ascension of the ascending node (RAAN) (rad)    [nx1]
%   w       -   Argument of perigee (rad)                             [nx1]
%   M       -   Mean anomaly (rad)                                    [nx1]
%
% =========================================================================
%
% Output:
%
%   CART    -   Position and velocity vector in ECI coordinates (km,  [nx6]
%               km/s)
%
% =========================================================================
%
% Initial version: Mar 2013;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

    sz = size(KEP);
    
    if sz(2) ~= 6 || ~all(isa(KEP, 'double'))
        error('Kep2Cart:InvalidOrbitalElements', 'Invalid Keplerian elements');
    end
    
    a     = KEP(:,1);
    e     = KEP(:,2);
    i     = KEP(:,3);
    Omega = KEP(:,4);
    w     = KEP(:,5);
    M     = KEP(:,6);
    
    if any(abs(e - 1) <= 1E-5)
        error('Kep2Cart:ParabolicOrRectilinearOrbit', 'Parabolic or rectilinear orbit, cannot calculate Cartesian elements');
    end
    
    if any(e > 1)
        error('Kep2Cart:HyperbolicOrbit', 'Hyperbolic orbit, cannot calculate Cartesian elements');
    end
    
    if any(abs(e) <= 1E-5)
        warning('Kep2Cart:SpecialCaseOrbit', 'Special case orbit: circular orbit');
    end
    
    if any(abs(i) <= 1E-5)
        warning('Kep2Cart:SpecialCaseOrbit', 'Special case orbit: equatorial orbit');
    end

    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu = 3.986004418e5;

    % Semilatus rectum
    p = a .* (1 - e.^2);
    
    % Compute eccentric anomaly from mean anomaly
    EA = Mean2EccAnomaly(M,e);
    
    % Compute true anomaly from eccentric anomaly
    nu = 2 * atan(sqrt((1+e)./(1-e)).*tan(EA/2));
    
    % Position and Velocity in Perifocal Coordinate System (PQW)
    x  = (p.*cos(nu)) ./ (1 + e.*cos(nu));
    y  = (p.*sin(nu)) ./ (1 + e.*cos(nu));
    z  = zeros(sz(1),1);
    vx = -sqrt(mu./p) .* sin(nu);
    vy =  sqrt(mu./p) .* (e + cos(nu));
    vz =  z;
    
    a11 = cos(Omega).*cos(w)-sin(Omega).*sin(w).*cos(i);
    a12 = -1*cos(Omega).*sin(w)-sin(Omega).*cos(w).*cos(i);
    a13 = sin(Omega).*sin(i);
    a21 = sin(Omega).*cos(w)+cos(Omega).*sin(w).*cos(i);
    a22 = -1*sin(Omega).*sin(w)+cos(Omega).*cos(w).*cos(i);
    a23 = -1*cos(Omega).*sin(i);
    a31 = sin(w).*sin(i);
    a32 = cos(w).*sin(i);
    a33 = cos(i);
    
    col1 = dot([a11,a12,a13],[x,y,z],2);
    col2 = dot([a21,a22,a23],[x,y,z],2);
    col3 = dot([a31,a32,a33],[x,y,z],2);
    col4 = dot([a11,a12,a13],[vx,vy,vz],2);
    col5 = dot([a21,a22,a23],[vx,vy,vz],2);
    col6 = dot([a31,a32,a33],[vx,vy,vz],2);
    
    CART = [col1,col2,col3,col4,col5,col6];
    
    CART(abs(CART) <= 1E-5) = 0;
    
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
% D. Plakalovic  | 12-11-2014 | Added more detailed commentary
% E. White       | 06-13-2023 | Added warnings for special case orbits,
%                               added errors for parabolic orbits, added
%                               input validation for size and type,
%                               updated documentation, added compliant
%                               header and footer
% E. White       | 06-30-2023 | Updated documentation, added checks for
%                               rectilinear and hyperbolic orbits
% 
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

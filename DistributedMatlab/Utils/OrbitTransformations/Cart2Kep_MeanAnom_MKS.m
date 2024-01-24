function [a,e,i,Omega,w,anom] = Cart2Kep_MeanAnom_MKS(r,v)
% Cart2Kep_MeanAnom_MKS - Converts Cartesian elements to Keplerian elements
%
% Syntax: [a,e,i,Omega,w,anom] = Cart2Kep_MeanAnom_MKS(r,v);
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
%   Converts Cartesian elements to Keplerian elements. The state vector 
%   must be input in units of m and m/s.
%   The function does not handle any special case orbits such as 
%   equatorial circular (True Longitude of Perigee needed), circular
%   inclined (Argument of Latitude needed), or equatorial elliptical (True 
%   Longitude needed).
%
%   NOTE: This is a less general version of the Cart2Kep function,
%   optimized for use in the K2Balt task (see Pc_MC_Kep2body_parallel); for
%   general-use purposes, please use Cart2Kep. DO NOT increase the
%   robustness and/or error-checking of this function, as it is optimized
%   for speed ONLY.
%
% =========================================================================
%
% Input:
%
%   r       -   Position vector in ECI coordinates (m)       [1x3] or [3x1]
%   v       -   Velocity vector in ECI coordinates (m/s)     [1x3] or [3x1]
%
% =========================================================================
%
% Output:
%
%   a       -   Semi-major axis (m)
%   e       -   Eccentricity
%   i       -   Inclination (rad [0, 2pi])
%   Omega   -   Right ascension of the ascending node (RAAN) (rad [0, 2pi])
%   w       -   Argument of perigee (rad [0, 2pi])
%   anom    -   Anomaly angle (mean) (rad [0, 2pi])  
%
% =========================================================================
%
% Initial version: Mar 2013;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

    % Earth gravitational constant (EGM-96) [m^3/s^2]
    mu = 3.986004418e14;

    % Angular momentum vector
    h = cross(r,v);
    
    % Line of nodes vector
    n = cross([0;0;1],h);
            
    % Magnitude of position, velocity, angular momentum, and line of nodes vectors
    norm_r = norm(r);
    norm_v = norm(v);
    norm_h = norm(h);
    norm_n = norm(n);
    
    % Specific energy
    Energy = 0.5 * norm_v^2 - mu/norm_r;
    
    % Semi-major axis
    a = -0.5 * mu / Energy;
    
    % Eccentricity: Vector and its norm
    E = cross(v,h)/mu - r/norm_r;
    e = norm(E);
      
    % Inclination  
    i = real(acos(h(3)/norm_h));
    
    % Right ascension of the ascending node (RAAN)
    Omega = real(acos(n(1)/norm_n));
    
    % RAAN quadrant check
    if (n(2) < 0)
        Omega = 2*pi - Omega;
    end
    
    % Argument of perigee
    w = real(acos(dot(n,E)/(norm_n*e)));
    
    % Argument of perigee quadrant check
    if (E(3) < 0)
        w = 2*pi - w;
    end
    
    % True anomaly
    anom = real(acos(dot(E,r)/(norm_r*e)));

    % True anomaly quadrant check
    if (dot(r,v) < 0)
        anom = 2*pi - anom;
    end

    % Eccentric anomaly
    anom  = 2 * atan(sqrt((1-e)/(1+e))*tan(anom/2));

    % Eccentric anomaly quadrant check
    if (anom < 0)
        anom = anom + 2*pi;
    end

    % Mean anomaly
    anom = (anom - e * sin(anom));
    
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
% D. Plakalovic  | 12-11-2014 |  Added functionality to include mean and 
%                                eccentric anomaly as an output option as
%                                well as an option to choose degrees or 
%                                radians as units of output
% D. Hall        | 2016-06-28 |  Revised to use MKS units, accept (r,v)
%                                column vectors, and return separated
%                                Keplerian elements.
% E. White       | 2023-06-12 |  Added compliant header and footer and
%                                updated documentation
% E. White       | 2023-06-13 |  Fixed error in documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

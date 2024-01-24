function [r,v] = Kep2Cart_MeanAnom_MKS(a,e,i,Omega,w,M)
% Kep2Cart_MeanAnom_MKS - Converts Keplerian elements to Cartesian elements
%
% Syntax: [r,v] = Kep2Cart_MeanAnom_MKS(a,e,i,Omega,w,M);
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
%   Converts Keplerian elements to Cartesian elements. Note that the input 
%   anomaly angle must be with respect to mean anomaly and all angle 
%   inputs are in the units of radians.
%
%   NOTE: Mean2EccAnomaly.m is a dependency for this function.
%
%   NOTE: This is a less general version of the Kep2Cart function,
%   optimized for use in the K2Balt task (see Pc_MC_Kep2body_parallel); for
%   general-use purposes, please use Kep2Cart. DO NOT increase the
%   robustness and/or error-checking of this function, as it is optimized
%   for speed ONLY.
%
% =========================================================================
%
% Input:
%
%   a       -   Semi-major axis (m)
%   e       -   Eccentricity
%   i       -   Inclination (rad)
%   Omega   -   Right ascension of the ascending node (RAAN) (rad)
%   w       -   Argument of perigee (rad)
%   M       -   Mean anomaly (rad)
%
% =========================================================================
%
% Output:
%
%   r       -   Position vector in ECI coordinates (m)                [3x1] 
%   v       -   Velocity vector in ECI coordinates (m/s)              [3x1] 
%
% =========================================================================
%
% Initial version: Mar 2013;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

    % Earth gravitational constant (EGM-96) [m^3/s^2]
    mu = 3.986004418e14;

    % Semilatus rectum
    p = a .* (1 - e.^2);
    
    % Compute eccentric anomaly from mean anomaly
    EA = Mean2EccAnomaly(M,e);
    
    % Compute true anomaly from eccentric anomaly
    nu = 2 * atan(sqrt((1+e)./(1-e)).*tan(EA/2));
    
    % Position and Velocity in Perifocal Coordinate System (PQW)
    
    cos_nu = cos(nu); sin_nu = sin(nu);
    ecos_nu = e.*cos_nu;
    sqrt_muop = sqrt(mu./p);
    
    x  = (p.*cos_nu) ./ (1 + ecos_nu);
    y  = (p.*sin_nu) ./ (1 + ecos_nu);    
    z  = zeros(size(x));
    
    vx = -sqrt_muop .* sin_nu;
    vy =  sqrt_muop .* (e + cos_nu);
    vz =  z;
    
    cos_O = cos(Omega); sin_O = sin(Omega);
    cos_w = cos(w); sin_w = sin(w);
    cos_i = cos(i); sin_i = sin(i);
    
    sin_w_cos_i = sin_w.*cos_i;
    cos_w_cos_i = cos_w.*cos_i;

    a11 = cos_O.*cos_w-sin_O.*sin_w_cos_i;
    a12 = -cos_O.*sin_w-sin_O.*cos_w_cos_i;
    a13 = sin_O.*sin_i;
    a21 = sin_O.*cos_w+cos_O.*sin_w_cos_i;
    a22 = -sin_O.*sin_w+cos_O.*cos_w_cos_i;
    a23 = -cos_O.*sin_i;
    a31 = sin_w.*sin_i;
    a32 = cos_w.*sin_i;
    a33 = cos_i;
    
    col1 = dot([a11,a12,a13],[x,y,z],2);
    col2 = dot([a21,a22,a23],[x,y,z],2);
    col3 = dot([a31,a32,a33],[x,y,z],2);
    col4 = dot([a11,a12,a13],[vx,vy,vz],2);
    col5 = dot([a21,a22,a23],[vx,vy,vz],2);
    col6 = dot([a31,a32,a33],[vx,vy,vz],2);
    
    r = [col1;col2;col3];
    v = [col4;col5;col6];
        
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
% D. Hall        | 2016-06-28 | Revised to use MKS units, accept separated
%                               Keplerian elements, and return (r,v)
%                               column vectors.
% E. White       | 2023-06-30 | Added compliant header and footer and
%                               updated documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

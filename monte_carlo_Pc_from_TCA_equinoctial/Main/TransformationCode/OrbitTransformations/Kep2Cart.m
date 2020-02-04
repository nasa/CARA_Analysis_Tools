function [CART] = Kep2Cart(KEP)
% function [r,v] = Kep2Cart(a,e,i,Omega,w,M)
%
% Kep2Cart - Converts an object's Keplerian elements to Cartesian elements.
%            Note that the input anomaly angle must be with respect to mean
%            anomaly, and all angles inputs are in the units of radians.
%            It does not handle any special case orbits such as 
%            equatorial circular, circular inclined, or equatorial elliptical.
%
% Syntax:    [r,v] = Kep2Cart(a,e,i,Omega,w,M)
%
% Inputs: 
%    a     - Semi-Major axis                               Units: km
%    e     - Eccentricity                                  Units: None
%    i     - Inclination                                   Units: rad   [-2pi,2pi]
%    Omega - Right Ascension of the Ascending Node (RAAN)  Units: rad   [-2pi,2pi]
%    w     - Argument of Perigee                           Units: rad   [-2pi,2pi]
%    M     - Mean Anomaly                                  Units: rad   [-2pi,2pi]
%
% Outputs:
%    r     - Position vector in ECI coordinates.           Units: km 
%            (1x3 row vector)
%    v     - Velocity vector in ECI coordinates.           Units: km/s
%            (1x3 row vector)
%
% Examples/Validation Cases:
%
% Other m-files required:  Mean2EccAnomaly.m
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% March 2013; Last Revision: 11-Dec-2014
%
% ----------------- BEGIN CODE -----------------

    % KEP = [a,e,i,Omega,w,M]

    sz = size(KEP,1);
    
    a     = KEP(:,1);
    e     = KEP(:,2);
    i     = KEP(:,3);
    Omega = KEP(:,4);
    w     = KEP(:,5);
    M     = KEP(:,6);

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
    z  = zeros(sz,1);
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
    
    
%     % Transformation Matrix from Perifocal to ECI (3 Euler rotations) 
%     PQW2ECI = [cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i), ...
%                -1*cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i), ...
%                sin(Omega)*sin(i); ...
%                sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i), ...
%                -1*sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i), ...
%                -1*cos(Omega)*sin(i); ...
%                sin(w)*sin(i), cos(w)*sin(i), cos(i)];
%     
%     % Convert position and velocity from PQW to ECI
%     r = transpose(PQW2ECI * [x,y,z]');
%     v = transpose(PQW2ECI * [vx,vy,vz]');
    
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
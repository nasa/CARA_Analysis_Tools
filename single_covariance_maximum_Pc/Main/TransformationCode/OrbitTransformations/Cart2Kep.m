function [KEP] = Cart2Kep(CART,AnomType,UnitType)
%
% Cart2Kep - Converts an object's Cartesian elements to Keplerian elements.
%            The state vector must be input in units of km and km/s. The 
%            user specifies which anomaly angle to be output (true, mean, 
%            eccentric) as well as the units of the output angles (degrees
%            or radians).
%            The function does not handle any special case orbits such as 
%            equatorial circular (True Longitude of Perigee needed), circular
%            inclined (Argument of Latitude needed), or equatorial elliptical 
%            (True Longitude needed).
%
% Syntax:        [a,e,i,Omega,w,anom] = Cart2Kep(r,v,'AnomType','UnitType')
%            Ex: [a,e,i,Omega,w,M]    = Cart2Kep(r,v,'Mean','Rad')
%
% Inputs: 
%    r     - Position vector in ECI coordinates.           Units: km 
%            (1x3 row vector or 3x1 column vector)
%    v     - Velocity vector in ECI coordinates.           Units: km/s
%            (1x3 row vector or 3x1 column vector)
% AnomType - Anomaly angle (String = 'True', 'Mean', or 'Eccentric')
% UnitType - Angles output units (String = 'Deg' or 'Rad')
%
% Outputs: 
%    a     - Semi-Major axis                               Units: km
%    e     - Eccentricity                                  Units: None
%    i     - Inclination                                   Units: deg or rad  [0,360] or [0,2pi]
%    Omega - Right Ascension of the Ascending Node (RAAN)  Units: deg or rad  [0,180] or [0, pi]
%    w     - Argument of Perigee                           Units: deg or rad  [0,360] or [0,2pi]
%    anom  - Anomaly angle (True, Mean, or Eccentric)      Units: deg or rad  [0,360] or [0,2pi]
%
% Examples/Validation Cases:
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% Revision History:
%   2/20/2018:            Added Accomodation of "0" state elements

% Latest Revision:        7/17/2019
%
% ----------------- BEGIN CODE -----------------

    % Ensure that cartesian states are non-zero
    for i=1:length(CART)
        if CART(i)==0
            CART(i)=1E-5;
        end
    end
    r = CART(:,1:3);
    v = CART(:,4:6);

    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu = 3.986004418e5;

    % Angular momentum vector
    h = cross(r,v);
    
    % Unit vector in k-direction
    k = [0,0,1];
    
    % Line of nodes vector
    n = cross(k,h);
            
    % Magnitude of position, velocity, angular momentum, and line of nodes vectors
    norm_r = norm(r); norm_v = norm(v); norm_h = norm(h); norm_n = norm(n);
    
    % Specific energy
    Energy = 0.5 * norm_v^2 - mu/norm_r;
    
    % Semi-major axis
    a = -0.5 * mu / Energy;
    
    % Eccentricity: Vector and its norm
    E = cross(v,h)/mu - r/norm_r;
    e = norm(E);
      
    % Inclination  
    i = real(acosd(h(3)/norm_h));
    
    % Right ascension of the ascending node (RAAN)
    Omega = real(acosd(n(1)/norm_n));
    
    % RAAN quadrant check
    if (n(2) < 0)
        Omega = 360 - Omega;
    end
    
    % Argument of perigee
    w = real(acosd(dot(n,E)/(norm_n*e)));
    
    % Argument of perigee quadrant check
    if (E(3) < 0)
        w = 360 - w;
    end
    
    % Compute anomaly angle
    switch lower(AnomType)
        
        case 'true'
    
            % True anomaly
            anom = real(acosd(dot(E,r)/(norm_r*e)));
    
            % True anomaly quadrant check
            if (dot(r,v) < 0)
                anom = 360 - anom;
            end
        
        case 'eccentric'
            
            % True anomaly
            anom = real(acos(dot(E,r)/(norm_r*e)));
    
            % True anomaly quadrant check
            if (dot(r,v) < 0)
                anom = 2*pi - anom;
            end
            
            % Eccentric anomaly
            anom  = (2 * atan(sqrt((1-e)/(1+e))*tan(anom/2))) * 180/pi;
            
            % Eccentric anomaly quadrant check
            if (anom < 0)
                anom = anom + 360;
            end
            
        case 'mean'
            
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
            anom = (anom - e * sin(anom)) * 180/pi;
            
        otherwise
            error([AnomType ' as AnomType is not supported...']);
                
    end
    
    % Convert output to specified units
    switch lower(UnitType)
        
        case 'deg'
          
            % Do nothing, already in degrees 
            
        case 'rad'
            
            % Convert output to radians
            i     = i     * pi/180; 
            Omega = Omega * pi/180;
            w     = w     * pi/180;
            anom  = anom  * pi/180;
        
        otherwise
            error([UnitType ' as UnitType is not supported...']);
            
    end
    
    KEP = [a,e,i,Omega,w,anom];
    
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
% T. Lechtenberg | 07-17-2019 |  Added more flexibility to options inputs
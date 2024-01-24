function [KEP] = Cart2Kep(CART,AnomType,UnitType)
% Cart2Kep - Converts an object's Cartesian elements to Keplerian elements
%
% Syntax: [a,e,i,Omega,w,anom] = Cart2Kep(CART,'AnomType','UnitType');
% Example: [a,e,i,Omega,w,M] = Cart2Kep([r_vec,v_vec],'Mean','Rad');
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
%   Converts an object's Cartesian elements to Keplerian elements. The 
%   state vector must be input in units of km and km/s. The user specifies
%   which anomaly angle to be output (true, mean, eccentric) as well as the
%   units of the output angles (degrees or radians).
%   The function does not handle any special case orbits such as 
%   equatorial circular (True Longitude of Perigee needed), circular
%   inclined (Argument of Latitude needed), or equatorial elliptical (True 
%   Longitude needed).
%
% =========================================================================
%
% Input:
%
%   CART        -   Cartesian state vector                   [1x6] or [6x1]
%   AnomType    -   Anomaly angle ('True', 'Mean', or 'Eccentric')
%   UnitType    -   Angle output units ('Deg' or 'Rad')
%
% =========================================================================
%
% Output:
%
%   a       -   Semi-major axis (km)
%   e       -   Eccentricity
%   i       -   Inclination (deg [0, 180] or rad [0, pi])
%   Omega   -   Right ascension of the ascending node (RAAN) (deg [0, 360]
%               or rad [0, 2pi])
%   w       -   Argument of perigee (deg [0, 360] or rad [0, 2pi])
%   anom    -   Anomaly angle (true, mean, or eccentric) (deg [0, 360] or 
%               rad [0, 2pi])
%
% =========================================================================
%
% Initial version: Mar 2013;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

    % Ensure that Cartesian elements are valid
    if ~all(isa(CART, 'double')) || length(CART) ~= 6
        error('Cart2Kep:InvalidCartesianElements', 'Input is not a valid set of Cartesian elements');
    end

    r = CART(1:3);
    v = CART(4:6);
    
    if all(~v)
        error('Cart2Kep:ZeroVelocity', 'Satellite cannot have zero velocity');
    end
    
    if norm(r) <= 6378.1370
        error('Cart2Kep:InvalidPosition', 'Satellite cannot be inside the Earth');
    end

    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu = 3.986004418e5;

    % Angular momentum vector
    h = cross(r,v);
    
    if abs(h) <= 1E-5
        error('Cart2Kep:RectilinearOrbit', 'Rectilinear orbit, unable to calculate Keplerian elements');
    end
    
    % Unit vector in k-direction
    k = [0,0,1];
    
    % Line of nodes vector
    n = cross(k,h);
            
    % Magnitude of position, velocity, angular momentum, and line of nodes vectors
    norm_r = norm(r); norm_v = norm(v); norm_h = norm(h); norm_n = norm(n);
    
    % Eccentricity: Vector and its norm
    E = cross(v,h)/mu - r/norm_r;
    e = norm(E);
    
    if abs(1 - e) <= 1E-5
        error('Cart2Kep:ParabolicOrbit', 'Parabolic orbit, unable to calculate Keplerian elements');
    end
    
    if e > 1
        error('Cart2Kep:HyperbolicOrbit', 'Hyperbolic orbit, unable to calculate Keplerian elements');
    end

    % Specific energy
    Energy = 0.5 * norm_v^2 - mu/norm_r;

    % Semi-major axis
    a = -0.5 * mu / Energy; 
      
    % Inclination  
    i = real(acosd(h(3)/norm_h));
    
    if abs(e) <= 1E-5 && abs(norm_n) <= 1E-5
        warning('Cart2Kep:EquatorialCircularOrbit', 'Special case orbit: equatorial circular orbit');
        
        e = 0;
        
        % No ascending node
        Omega = 0;
        
        % No perigee
        w = 0;
        
        % True longitude
        anom = real(acosd(r(1)/norm_r));
        
        % True longitude quadrant check
        if r(2) < 0
            anom = 360 - anom;
        end
    elseif abs(e) <= 1E-5
        warning('Cart2Kep:InclinedCircularOrbit', 'Special case orbit: inclined circular orbit');
        
        e = 0;
        
        % Right ascension of the ascending node (RAAN)
        Omega = real(acosd(n(1)/norm_n));

        % RAAN quadrant check
        if (n(2) < 0)
            Omega = 360 - Omega;
        end

        % No perigee
        w = 0;
        
        % Argument of latitude
        anom = real(acosd(dot(n,r)/(norm_n*norm_r)));
        
        % Argument of latitude unit check
        if r(3) < 0
            anom = 360 - anom;
        end
    elseif abs(norm_n) <= 1E-5
        warning('Cart2Kep:GeneralEquatorialOrbit', 'Special case orbit: equatorial elliptical or hyperbolic orbit');

        % No ascending node
        Omega = 0;

        % True longitude of periapsis
        w = real(acosd(E(1)/e));

        % True longitude of periapsis quadrant check
        if (E(2) < 0)
            w = 360 - w;
        end
        
        % True anomaly
        anom = real(acosd(dot(E,r)/(norm_r*e)));

        % True anomaly quadrant check
        if (dot(r,v) < 0)
            anom = 360 - anom;
        end
    else
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

        % True anomaly
        anom = real(acosd(dot(E,r)/(norm_r*e)));

        % True anomaly quadrant check
        if (dot(r,v) < 0)
            anom = 360 - anom;
        end
    end
    
    % Compute anomaly angle
    switch lower(AnomType)
        
        case 'true'
    
            % Do nothing, already true
        
        case 'eccentric'
            
            % Eccentric anomaly
            anom  = (2 * atand(sqrt((1-e)/(1+e))*tand(anom/2)));
            
            % Eccentric anomaly quadrant check
            if (anom < 0)
                anom = anom + 360;
            end
            
        case 'mean'
            
            % Eccentric anomaly
            anom  = 2 * atand(sqrt((1-e)/(1+e))*tand(anom/2));
            
            % Eccentric anomaly quadrant check
            if (anom < 0)
                anom = anom + 360;
            end
            
            % Mean anomaly
            
            % Since eccentricity is unitless, the calculation with sin ends
            % up being in radians, so the rad2deg adjustment is necessary
            % here
            anom = (anom - rad2deg(e * sind(anom)));
            
        otherwise
            error('Cart2Kep:InvalidAnomType', [AnomType ' as AnomType is not supported...']);
                
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
            error('Cart2Kep:InvalidUnitType', [UnitType ' as UnitType is not supported...']);
            
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
% D. Plakalovic  | Mar - 2013 | Initial Development
% D. Plakalovic  | 12-11-2014 | Added functionality to include mean and 
%                               eccentric anomaly as an output option as
%                               well as an option to choose degrees or 
%                               radians as units of output
% T. Lechtenberg | 07-17-2019 | Added more flexibility to options inputs
% E. White       | 06-12-2023 | Added support for special case orbits,
%                               added check for parabolic orbits, added
%                               compliant header and footer, added input
%                               verification
% E. White       | 06-30-2023 | Fixed bug involving negative positional
%                               values, updated documentation, added
%                               checks for rectilinear and hyperbolic
%                               orbits

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

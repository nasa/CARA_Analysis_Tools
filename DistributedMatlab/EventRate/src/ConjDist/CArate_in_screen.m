function [CA_in_screen, LT_in_screen] = CArate_in_screen( ...
    Dbin, DbinSquared, in_bin, rspRIC, vspRIC, params)
% CArate_in_screen - Calculates which spherical screening volumes are
%                    penetrated by a conjuncting secondary
% Syntax: [CA_in_screen, LT_in_screen] = ...
%    CArate_in_screen(Dbin, DbinSquared, in_bin, rspRIC, vspRIC, params)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculate which inscribed screening volumes have been 
% penetrated during a CA event.  These screening volumes are inscribed into 
% the set of spheres defined by the CA distance bins.
%
% This function uses a linear approximation for the trajectory of the
% secondary in the primary's RIC frame.
%
% =========================================================================
%
% Input:
%
%   Dbin        - The radii of the spherical bins (km).               [1xN]
%
%   DbinSquared - Dbin.^2                                             [1xN]
%
%   in_bin      - A logical array indicating which spherical bins     [1xN]
%                 contain the CA point.
%
%   rspRIC      - The CA position vector of the secondary in the      [3x1]
%                 primary's RIC reference frame.
%
%   rspRIC      - The CA velocity vector of the secondary in the      [3x1]
%                 primary's RIC reference frame.
%
%   params      - CArate parameters; see function 
%                 "CArate_default_params."
%
% =========================================================================
%
% Output:
%
%   CA_in_screen - A logical array indicating which inscribed         [1xN]
%                  non-spherical screening volumes contain the 
%                  CA point.
%
%   LT_in_screen - A logical array indicating which inscribed         [1xN]
%                  non-spherical screening volumes have been 
%                  penetrated by the linearized trajectory.
%
%
% =========================================================================
%
% Dependencies:
%
%   None
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

% Initialize output to indicate no screening volume incursions

CA_in_screen = false(size(in_bin));
LT_in_screen = CA_in_screen;

% Number of counts in the spherical bins

sum_in_bin = sum(in_bin);

% Incursions into the screening volumes can only occur if there are
% some incursions into the circumscribing spheres.

if sum_in_bin == 0
    return;
end

% Process different screening volume types

switch params.screen_type
    
    case 'Sphere'
        
        % Spherical screening volumes are simple because they coincide with
        % the circumscribing spherical bins
        
        CA_in_screen(in_bin) = true;
        LT_in_screen(in_bin) = true;

    case 'RIC_ellipsoid'
        
        % Check if secondary is within primary-centered screening ellipsoid
        
        find_in_bin = find(in_bin);

        Dbin2 = DbinSquared(in_bin);
        rspRICt = rspRIC';
        cc0 = params.Cnorm * rspRIC;
        cc0 = rspRICt * cc0;

        CvspRIC = params.Cnorm * vspRIC;
        aa = vspRIC' * CvspRIC;
        bbhalfsq = (rspRICt * CvspRIC)^2;

        for nnbin=1:sum_in_bin

            % Current bin number

            nbin = find_in_bin(nnbin);

            if (cc0 <= Dbin2(nnbin))
                % CA point itself is inside the ellipsoid.
                CA_in_screen(nbin) = true;
                LT_in_screen(nbin) = true;
            else
                % Check if the linearized trajectory penetrates
                % the ellipsoid.  This can be done using a quadratic
                % equation analysis where the descriminant is nonnegative.
                cc = cc0 - Dbin2(nnbin);
                if (bbhalfsq >= aa*cc)
                    LT_in_screen(nbin) = true;
                end
            end

        end
        
    case 'RIC_box'

        % Check if secondary is within primary-centered screening box
        
        find_in_bin = find(in_bin);

        for nnbin=1:sum_in_bin
            
            % Current bin number
            
            nbin = find_in_bin(nnbin);
            
            % Half side-lengths of the box inscribed into the bin
            % sphere
            
            abc = Dbin(nbin) * params.screen_ratios;
            
            % Check if the CA point itself is inside the box,
            % otherwise check if the linearized trajectory intersects
            % the surface.
            
            if (abs(rspRIC(1)) <= abc(1)) && ...
               (abs(rspRIC(2)) <= abc(2)) && ...
               (abs(rspRIC(3)) <= abc(3))
           
                % CA point itself is inside the box
                
                CA_in_screen(nbin) = true;
                LT_in_screen(nbin) = true;
                
            else
                
                % Check for intersections of the linearized trajectory
                % with any of the six box faces
                
                for nf=1:6
                    
                    switch nf
                        case 1
                            % +x face
                            i=1; j=2; k=3;
                            pm = +1;
                        case 2
                            % -x face
                            pm = -1;
                        case 3
                            % +y face
                            i=2; j=1; k=3;
                            pm = +1;
                        case 4
                            % -y face
                            pm = -1;
                        case 5
                            % +z face
                            i=3; j=2; k=1;
                            pm = +1;
                        case 6
                            % -z face
                            pm = -1;
                    end
                    
                    % Transit time along linearized trajectory to where it
                    % intersects the plane corresponding to the current
                    % box face
                    
                    tau = (pm*abc(i)-rspRIC(i))/vspRIC(i);
                    
                    % Intersection point
                    
                    rtau = rspRIC+vspRIC*tau;
                    
                    % If the intersection point is within the face
                    % boundaries, then a box incursion has occurred,
                    % and there is no need to check for others
                    
                    if (abs(rtau(j)) <= abc(j)) && (abs(rtau(k)) <= abc(k))
                        LT_in_screen(nbin) = true;
                        break;
                    end
                    
                end

            end
            
        end

    otherwise

        error(['Unrecognized or unimplemented screening volume ' ...
               'type: ' params.screen_type]);

end

return;
end
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
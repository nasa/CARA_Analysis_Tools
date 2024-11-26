function lat = SameSineLatitude(lat,InRadians)
% SameSineLatitude - Return the angle within latitude limits with the same 
%                    sine value.
% Syntax: lat = SameSineLatitude(lat,InRadians);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

if InRadians
    halfpi = pi/2;
    if lat < -halfpi
        lat = -pi-lat;
    elseif lat > halfpi
        lat =  pi-lat;
    end    
else
    if lat < -90
        lat = -180-lat;
    elseif lat > 90
        lat =  180-lat;
    end
end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall    |  2024-Nov  | Initial version.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
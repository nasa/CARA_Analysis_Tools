function [x, y] = GenerateEllipsePoints(a, b, angDeg, shift)
% GenerateEllipsePoints - Generates 2D x-y coordinates to plot an ellipse
%                         according to the parameters passed in.
%
% Syntax: [x, y] = GenerateEllipsePoints(a, b, angDeg, shift);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Generates a set of x-y coordinates in order to graph an ellipse on a 2D
%   plane. Input parameters include the semi-major axis, semi-minor axis,
%   rotation (counter-clockwise) of the ellipse, and a translation of the
%   center of the ellipse.
%
%   The points are generated in such a way that the vertices are always
%   included as an x-y coordinate.
%
% =========================================================================
%
% Input:
%
%    a - (Optional) Semi-major axis of the ellipse                    [1x1]
%        Default is 1
%    b - (Optional) Semi-minor axis of the ellipse                    [1x1]
%        Default is the value of a
%    angDeg - (Optional) Clockwise angular rotation of the ellipse.   [1x1]
%             A 0 degree angDeg aligns the semi-major axis with the
%             x-axis.
%             Default is 0
%    shift - (Optional) x-y coordinates of the center of the ellipse. [1x2]
%            Default is [0 0]
%
% =========================================================================
%
% Output:
%
%    x - x-axis values of the perimeter of the ellipse             [1x1441]
%    y - y-axis values of the perimeter of the ellipse             [1x1441]
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Mar 2023
%
% ----------------- BEGIN CODE -----------------

    % Set default values based on optional inputs
    if nargin == 0
        a = 1;
        b = a;
        angDeg = 0;
        shift = [0 0];
    elseif nargin == 1
        b = a;
        angDeg = 0;
        shift = [0 0];
    elseif nargin == 2
        angDeg = 0;
        shift = [0 0];
    elseif nargin == 3
        shift = [0 0];
    elseif nargin ~= 4
        error('Invalid number of arguments passed in');
    end
    
    % Check the dimensions of the input parameters
    if ~isequal(size(a),[1 1])
        error('Semi-major axis must be a scalar value');
    end
    if ~isequal(size(b),[1 1])
        error('Semi-minor axis must be a scalar value');
    end
    if ~isequal(size(angDeg),[1 1])
        error('Angle (deg) must be a scalar value');
    end
    if ~isequal(size(shift),[1 2])
        error('Shift must be a 1x2 vector');
    end
    
    % Generate a nominal ellipse centered on the origin with the semi-major
    % axis aligned along the x-axis
    u = linspace(0,360,1441) * pi/180;
    x0 = a * cos(u);
    y0 = b * sin(u);
    
    % Rotate the ellipse by angDeg
    ang = angDeg * pi/180;
    x = x0.*cos(ang) - y0.*sin(ang);
    y = x0.*sin(ang) + y0.*cos(ang);
    
    % Move the ellipse to be centered on the shift coordinates
    x = x + shift(1);
    y = y + shift(2);
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-20 | Initial Development
% L. Baars       | 2025-Aug-25 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

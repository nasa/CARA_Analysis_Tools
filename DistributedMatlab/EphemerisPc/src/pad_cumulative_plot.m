function [x,y] = pad_cumulative_plot(xrng,x,y)
% pad_cumulative_plot - Pad a cumulative plot.
%
% Syntax: [x, y] = pad_cumulative_plot(xrng, x, y);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    xrng       -   X-axis range [xmin xmax]                [1x2]
%
%    x          -   X data array                            [1xN]
%
%    y          -   Y data array                            [1xN]
%
% =========================================================================
%
% Output:
%
%   x           -   Padded X data array
%
%   y           -   Padded Y data array
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

n = numel(x);
x = reshape(x,[1 n]);
y = reshape(y,[1 n]);

if (x(1) > xrng(1))
    x = cat(2,xrng(1),x);
    y = cat(2,0,y);
end
if (x(end) < xrng(2))
    x = cat(2,x,xrng(2));
    y = cat(2,y,y(end));
end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
function PD2 = interp_PD2(t,eph1,eph2,NoWeighting)
% interp_PD2 - Interpolate the squared relative positional distance (PD2)
%              between two ephemeris objects.
%
% Syntax: PD2 = interp_PD2(t, eph1, eph2);
%         PD2 = interp_PD2(t, eph1, eph2, NoWeighting);
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
%    t          -   Times to interpolate             [1xM] or [Mx1]
%
%    eph1       -   Primary ephemeris structure with fields:
%                     .T - Ephemeris times            [1xN]
%                     .X - State vectors              [6xN]
%                     .N - Number of points
%
%    eph2       -   Secondary ephemeris structure (same fields as eph1)
%
%    NoWeighting -  Flag to suppress weighted interpolation
%                   (optional, default = false)
%
% =========================================================================
%
% Output:
%
%   PD2         -   Squared relative positional distance   [1xM] or [Mx1]
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Default is to use weighted interpolation for smooth interpolant curves
if nargin < 4; NoWeighting = []; end
if isempty(NoWeighting); NoWeighting = false; end

% Allocate the output array
PD2 = NaN(size(t));

% Calculate the relative distances
M = numel(t);
for m=1:M
    % Interpolate both ephemeris tables
    r1 = interpState(t(m),eph1,NoWeighting);
    r2 = interpState(t(m),eph2,NoWeighting);
    % Calculate the RD value
    PD2(m) = sum((r2-r1).^2);
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

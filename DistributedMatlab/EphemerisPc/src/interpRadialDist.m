function rdist = interpRadialDist(t,eph,NoWeighting)
% interpRadialDist - Interpolation of radial distance magnitudes.
%
% Syntax: rdist = interpRadialDist(t, eph);
%         rdist = interpRadialDist(t, eph, NoWeighting);
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
%    eph        -   Ephemeris structure with fields:
%                     .T - Ephemeris times            [1xN]
%                     .X - State vectors              [6xN]
%                     .N - Number of points
%
%    NoWeighting -  Flag to suppress weighted interpolation
%                   (optional, default = false)
%
% =========================================================================
%
% Output:
%
%   rdist       -   Radial distance magnitudes       [1xM] or [Mx1]
%                                                    (same as input t)
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Default is to use weighted interpolation for smooth interpolant curves
if nargin < 3; NoWeighting = []; end
if isempty(NoWeighting); NoWeighting = false; end

% Calculate radial distances
rdist = NaN(size(t));
N = numel(t);
for n=1:N
    rdist(n) = sqrt(sum(interpState(t(n),eph,NoWeighting).^2));
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
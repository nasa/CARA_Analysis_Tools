function [c,b,coincidence] = find_bracketing_indices(T,Teph)
% find_bracketing_indices - Find the closest ephemeris point and
%                           bracketing ephemeris point, assuming that
%                           the time is within the eph. bounds, i.e.,
%                           Teph(1) <= T <= Teph(end).
%
% Syntax: [c, b, coincidence] = find_bracketing_indices(T, Teph);
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
%    T          -   Time to bracket (must satisfy
%                   Teph(1) <= T <= Teph(end))
%
%    Teph       -   Array of ephemeris times        [1xN] or [Nx1]
%
% =========================================================================
%
% Output:
%
%   c           -   Index of closest ephemeris time
%
%   b           -   Index of bracketing ephemeris time
%
%   coincidence -   Flag indicating T exactly coincides with
%                   an ephemeris point
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

Neph = numel(Teph);

if (T < Teph(1)) || (T > Teph(Neph))
    error('Time not within ephemeris bounds');
end

% Find closest eph time index, c
[~,c] = min(abs(T-Teph));

% Find bracketing eph time index, b
if T > Teph(c)
    % Bracketing point is later
    b = c+1; coincidence = false;
elseif T < Teph(c)
    % Bracketing point is earlier
    b = c-1; coincidence = false;
else
    % Bracketing time coincidence with closest time
    coincidence = true;
    % b = c;
    % Pick the closest of the adjacent times as the bracketing point
    if c == 1
        b = 2;
    elseif c == Neph
        b = Neph-1;
    else
        m = c-1; p = c+1;
        if abs(Teph(c)-Teph(m)) <= abs(Teph(c)-Teph(p))
            b = m;
        else
            b = p;
        end
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
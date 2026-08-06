function [S,M2,M2HBR] = calc_sep_maha_distances(EPH,HBR,Lclip)
% calc_sep_maha_distances - Calculate separation and (optionally)
%                           Mahalanobis distances from the provided
%                           ephemeris table.
%
% Syntax: [S, M2, M2HBR] = calc_sep_maha_distances(EPH, HBR, Lclip);
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
%    EPH        -   Ephemeris table structure with fields:
%                     .T   - Time vector             [1xN] or [Nx1]
%                     .X1  - Primary state vectors   [3xN] or [6xN]
%                     .X2  - Secondary state vectors [3xN] or [6xN]
%                     .P1  - Primary covariances     [3x3xN] or [6x6xN]
%                     .P2  - Secondary covariances   [3x3xN] or [6x6xN]
%                     (Only the position portion of state vectors and 
% %                   covariance matrices are used)
%    HBR        -   Combined hard body radius
%    Lclip      -   Clipping limit for the eigenvalues
%
% =========================================================================
%
% Output:
%
%   S           -   Separation distance                     [1xN]
%
%   M2          -   Maha distance at center of coll sphere  [1xN]
%
%   M2HBR       -   Min Maha distance on coll sphere        [1xN]
%                   (approximated)
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

N = numel(EPH.T);

% Separation distance
r = EPH.X2(1:3,:)-EPH.X1(1:3,:);
S = sqrt(sum(r.*r,1));

if nargout < 2
    return
end

% Maha distances
M2 = NaN(size(S)); % Maha distance at center of coll sphere
M2HBR = M2;        % Min Maha distance on coll sphere (approximated)
for n=1:N
    A = EPH.P1(1:3,1:3,n)+EPH.P2(1:3,1:3,n);
    [M2(n),M2HBR(n)] = calc_maha_distances(r(:,n),A,HBR,Lclip);
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
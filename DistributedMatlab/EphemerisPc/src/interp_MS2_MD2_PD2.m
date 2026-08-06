function [MS2,MD2,PD2] = interp_MS2_MD2_PD2(T,eph1,eph2,HBR,Lclip)

% interp_MS2_MD2_PD2 - Interpolate square of relative position separation
%                       distance (PD2), the mean-state Maha. distance at
%                       the center of collision sphere (MD2), and
%                       approximated minimum Maha. distance on the
%                       collision sphere (MS2).
%
% Syntax: [MS2, MD2, PD2] = interp_MS2_MD2_PD2(T, eph1, eph2, HBR, Lclip);
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
%    T          -   Times to interpolate             [1xM] or [Mx1]
%
%    eph1       -   Primary ephemeris structure with fields:
%                     .T - Ephemeris times            [1xN]
%                     .X - State vectors              [6xN]
%                     .P - Covariance matrices        [6x6xN]
%                     .N - Number of points
%
%    eph2       -   Secondary ephemeris structure (same fields as eph1)
%
%    HBR        -   Combined hard body radius (km)
%
%    Lclip      -   Clipping limit for the eigenvalues
%
% =========================================================================
%
% Output:
%
%   MS2         -   Approximated minimum Maha. distance    [1xM] or [Mx1]
%                   on the collision sphere
%
%   MD2         -   Mean-state Maha. distance at the       [1xM] or [Mx1]
%                   center of collision sphere
%
%   PD2         -   Square of relative position            [1xM] or [Mx1]
%                   separation distance
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

M = numel(T);
S = size(T);
PD2 = NaN(S);
MD2 = NaN(S);
MS2 = NaN(S);

for m=1:M
    
    % Interpolate both ephemeris tables
    [r1,~,P1] = interpStateCov(T(m),eph1);
    [r2,~,P2] = interpStateCov(T(m),eph2);
        
    % Calculate the MD2HBR value
    r = (r2-r1)';
    A = P1(1:3,1:3)+P2(1:3,1:3);
    [MS2(m),MD2(m),PD2(m)] = calc_PD2_MD2_MS2(r,A,HBR,Lclip);
    
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
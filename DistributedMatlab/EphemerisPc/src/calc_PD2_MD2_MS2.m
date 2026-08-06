function [MS2,MD2,PD2] = calc_PD2_MD2_MS2(r,A,HBR,Lclip)
% calc_PD2_MD2_MS2 - Calculate square of relative position separation
%                    distance (PD2), the mean-state Maha. distance at the
%                    center of collision sphere (MD2), and approximated
%                    minimum Maha. distance on the collision sphere (MS2).
%
% Syntax: [MS2, MD2, PD2] = calc_PD2_MD2_MS2(r, A, HBR, Lclip);
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
%    r          -   Relative position vector                [3x1]
%    A          -   Relative position covariance matrix     [3x3] or [6x6]
%                   (only upper-left 3x3 block is used)
%    HBR        -   Combined hard body radius
%    Lclip      -   Clipping limit for the eigenvalues
%
% =========================================================================
%
% Output:
%
%   MS2         -   Approximated minimum Maha. distance on the collision 
%                   sphere
%
%   MD2         -   Mean-state Maha. distance at the center of collision 
%                   sphere
%
%   PD2         -   Square of relative position separation distance
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Relative positional distance squared
PD2 = r'*r;

% Invert the rel. pos. covariance using eigenvalue clipping algorithm to 
% handle the (rare) occurrences of NPD matrices
[~,~,~,~,~,~,Ainv] = CovRemEigValClip(A(1:3,1:3),Lclip);
Air = Ainv * r;

% MD2 for center of collision sphere w/ radius HBR; this should be
% nonnegative
MD2 = r' * Air;

% Approximate the minimum Maha. distance over collision sphere, 
% to 1st order accuracy in the quantity HBR/|r|.
% Note: For HBR > 0, this approximation can yield negative estimates
if HBR > 0
    MS2 = MD2 - 2*HBR*norm(Air);
else
    MS2 = MD2;
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
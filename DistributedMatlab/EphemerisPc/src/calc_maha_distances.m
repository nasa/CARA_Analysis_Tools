function [M2,M2HBR] = calc_maha_distances(r,A,HBR,Lclip)
% calc_maha_distances - Calculate Mahalanobis distances from relative
%                       position vector and relative position covariance
%                       matrix.
%
% Syntax: [M2, M2HBR] = calc_maha_distances(r, A, HBR, Lclip);
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
%   M2          -   Mahalanobis distance squared for the center of the
%                   collision sphere
%
%   M2HBR       -   Minimum squared Mahalanobis distance over the collision
%                   sphere
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Invert relative position 3x3 covariance matrix
% (remediating and NPD occurances using eigenvalue clipping method)
[~,~,~,~,~,~,Ainv] = CovRemEigValClip(A(1:3,1:3),Lclip);
Air = Ainv * r;

% MD2 for center of collision sphere w/ radius HBR;
% this should always be nonnegative
M2 = r' * Air;

% Min. MD2 over collision sphere, to 1st order accuracy in HBR/|r|;
% this approximation can yield negative estimates
M2HBR = M2 - 2*HBR*norm(Air);

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
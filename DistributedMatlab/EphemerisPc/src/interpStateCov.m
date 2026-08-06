function [ri,vi,Pi,na,ra,va,Pa,nb,rb,vb,Pb] = interpStateCov(T,eph,NoWeighting)
% interpStateCov - Interpolate a position/velocity state and associated
%                  covariance matrix from an ephemeris
%
% Syntax: [ri, vi, Pi, na, ra, va, Pa, nb, rb, vb, Pb] = interpStateCov(T, eph);
%         [ri, vi, Pi, na, ra, va, Pa, nb, rb, vb, Pb] = interpStateCov(T, eph, NoWeighting);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Interpolate a position/velocity state and associated covariance matrix 
%   from an ephemeris table at the single scalar time T.
%
%   By default, for points not near the edges of the ephemeris span,
%   this function performs "weighted" Lagrange interpolation.
%   Specifically, this means that it provides a weighted average of two
%   interpolations:
%      "a" uses the 5 points centered on the nearest eph. time, and
%      "b" uses the 5 points centered on the 2nd nearest eph. time.
%   Typically these two eph. times "bracket" the interpolated time,
%   satisfying t(na) <= T <= t(nb).
%   The default position (or velocity) interpolant is a weighted sum
%      ri = (1-w)*ra + w*rb   with   0 <= w <= 1   (similarly for vi
%      and Pi)
%   with ra denoting the 5-pt Lagrange interpolant centered on point a,
%   and  rb denoting the 5-pt Lagrange interpolant centered on point b.
%
%   This weighting algorithm provides continuous interpolant curves
%   r(T), v(T) and P(T). If no weighting is used, then these curves
%   have discontinuities near the mid-points of the eph. intervals, due
%   to the 5-pt set of points used in the Lagrange interpolation
%   changing discretely. The differences between the a and b
%   interpolants also provide an indication of the interpolation
%   uncertainty.
%
% =========================================================================
%
% Input:
%
%    T          -   Scalar time for the interpolation
%
%    eph        -   Ephemeris structure with fields:
%                     .T - Ephemeris times              [1xN]
%                     .X - State vectors                [6xN]
%                     .P - Covariance matrices          [6x6xN]
%                     .N - Number of points
%
%    NoWeighting -  Flag to suppress weighting algorithm
%                   (optional, default = false)
%
% =========================================================================
%
% Output:
%
%   (ri,vi,Pi)  -   Weighted pos/vel/cov interpolants for time T:
%                   [1x3], [1x3], and [6x6]
%
%   (na,nb)     -   Indices for the bracketing ephemeris points a and b.
%                   Point a corresponds to the closest ephemeris time to
%                   T. If T exactly coincides with an ephemeris point,
%                   T = t(na), then nb = na.
%
%   (ra,va,Pa)  -   Interpolated pos/vel/cov at time T, calculated
%                   using Lagrange interpolation with the 5-pt set
%                   {na-2 : na+2}: [1x3], [1x3], and [6x6]
%
%   (rb,vb,Pb)  -   Interpolated pos/vel/cov at time T, calculated
%                   using Lagrange interpolation with the 5-pt set
%                   {nb-2 : nb+2}: [1x3], [1x3], and [6x6]
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------


% Default is to use weighted interpolation for smooth interpolant curves
if nargin < 3; NoWeighting = []; end
if isempty(NoWeighting); NoWeighting = false; end

% Find best five bracketing eph points for interpolation
[~,na] = min(abs(T-eph.T));
n1 = na-2;
n2 = na+2;
if (n1 < 1)
    n1 = 1; n2 = 5;
elseif (n2 > eph.N)
    n1 = eph.N-4; n2 = eph.N;
end

% Interpolate state and covariance
i = n1:n2;
[ra,va,Pa] = StateCovInterp(T, ...
    eph.T(i), eph.X(1:3,i)', eph.X(4:6,i)', eph.P(:,:,i));

% If appropriate, interpolate states with the alternate five bracketing
% points, and create weighted average states
if NoWeighting || (T == eph.T(na))
    nb = na;
    rb = ra; vb = va; Pb = Pa;
    ri = ra; vi = va; Pi = Pa;
else
    % If not at eph point, generate the alternate five bracketing points
    if T > eph.T(na)
        nb = na+1;
    else
        nb = na-1;
    end
    k1 = nb-2;
    k2 = nb+2;
    if (k1 < 1)
        k1 = 1; k2 = 5;
    elseif (k2 > eph.N)
        k1 = eph.N-4; k2 = eph.N;
    end
    if (k1 == n1)
        % Alternate bracketing points are same as the original five
        nb = na;
        rb = ra; vb = va; Pb = Pa;
        ri = ra; vi = va; Pi = Pa;
    else
        % Alternate bracketing points are different than the original,
        % so perform the interpolation
        i = k1:k2;
        [rb,vb,Pb] = StateCovInterp(T, ...
            eph.T(i), eph.X(1:3,i)', eph.X(4:6,i)', eph.P(:,:,i));
        % Weight the two alternate sets of states by time differences
        w = (T-eph.T(na)) / (eph.T(nb)-eph.T(na));
        if (w <= 0) || (w >= 1)
            error('Invalid weight');
        end
        omw = 1-w;
        ri = omw*ra+w*rb; vi = omw*va+w*vb; Pi = omw*Pa+w*Pb;
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
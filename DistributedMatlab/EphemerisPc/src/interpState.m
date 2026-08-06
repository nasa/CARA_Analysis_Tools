function [ri,vi,na,ra,va,nb,rb,vb] = interpState(T,eph,NoWeighting,ExtrapPars)

% interpState - Interpolate a position/velocity state from an ephemeris.
%
% Syntax: [ri, vi] = interpState(T, eph);
%         [ri, vi] = interpState(T, eph, NoWeighting);
%         [ri, vi] = interpState(T, eph, NoWeighting, ExtrapPars);
%         [ri, vi, na, ra, va, nb, rb, vb] = interpState(__);
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
%   Interpolate a position/velocity state from an ephemeris table at the
%   single scalar time T.
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
%      ri = (1-w)*ra + w*rb   with   0 <= w <= 1   (similarly for vi)
%   with ra denoting the 5-pt Lagrange interpolant centered on point a,
%   and  rb denoting the 5-pt Lagrange interpolant centered on point b.
%
%   This weighting algorithm provides continuous interpolant curves
%   r(T) and v(T). If no weighting is used, then the interpolant curves
%   have discontinuities near the mid-points of the eph. intervals, due
%   to the 5-pt set of points used in the Lagrange interpolation
%   changing discretely. The differences between the a and b
%   interpolants also provide an indication of the interpolation
%   uncertainty.
%
%   This function also (optionally) allows state extrapolation, using
%   the approximation of constant velocity, rectilinear motion beyond
%   the eph. bounds. By default, this capability is not allowed, and is
%   not recommended except for short duration extrapolations.
%
% =========================================================================
%
% Input:
%
%    T          -   Scalar time for the interpolation
%
%    eph        -   Ephemeris structure with fields:
%                     .T             - Ephemeris times      [1xN]
%                     .X             - State vectors        [6xN]
%                     .N             - Number of points
%                     .SecPerTimeUnit - Seconds per time
%                                      unit (optional,
%                                      default = 86400)
%
%    NoWeighting -  Flag to suppress weighting algorithm
%                   (optional, default = false)
%
%    ExtrapPars -   Extrapolation parameters
%                   (optional, default = 0, no extrapolation)
%
% =========================================================================
%
% Output:
%
%   (ri,vi)     -   Weighted position/velocity state interpolants for
%                   time T, each [1x3]
%
%   (na,nb)     -   Indices for the bracketing ephemeris points a and b.
%                   Point a corresponds to the closest ephemeris time to
%                   T. If T exactly coincides with an ephemeris point,
%                   T = t(na), then nb = na.
%
%   (ra,va)     -   Interpolated position/velocity state at time T,
%                   calculated using Lagrange interpolation with the
%                   5-pt set {na-2 : na+2}, each [1x3]
%
%   (rb,vb)     -   Interpolated position/velocity state at time T,
%                   calculated using Lagrange interpolation with the
%                   5-pt set {nb-2 : nb+2}, each [1x3]
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Default is to use weighted interpolation for smooth interpolant curves
Nargin = nargin;
if Nargin < 3; NoWeighting = []; end
if isempty(NoWeighting); NoWeighting = false; end

% Default state extrapolation parameters
%  ExtrapPars - Allows state extrapolation outside of eph. time bounds:
%       ExtrapPars <= 0      =>  No extrapolation allowed; attempting
%                                to do so results in an error
%                                (default).
%       ExtrapPars = N > 0   =>  Extrapolation allowed within N
%                                timesteps of ephemeris endpoints with
%                                a warning; attempting longer
%                                extrapolations results in an error.
%       ExtrapPars = [M, N]  =>  Extrapolation allowed within M
%                                timesteps of allowed with no warning;
%                                extrapolation within M to N timesteps
%                                performed with a warning; 
%                                extrapolation beyond N timesteps
%                                results in an error.
%       Notes: Setting N = [Inf Inf] allows arbitrarily long extrapolations
%              with no warnings or errors.

if Nargin < 4; ExtrapPars = []; end
if isempty(ExtrapPars); ExtrapPars = 0; end
if all(ExtrapPars <= 0)
    ExtrapolationAllowed = false;
else
    ExtrapolationAllowed = true;
    Nex = numel(ExtrapPars);
    if Nex == 1
        % Issue warnings for any extrapolations
        ExtrapPars = [0 ExtrapPars];
    elseif Nex ~= 2
        error('Incorrect ExtrapPars parameter dimensions');
    % elseif ExtrapPars(1) >= ExtrapPars(2)
    %     warning('ExtrapPars(1) >= ExtrapPars(2) is not a sensible way to use the extrapolation mode');
    end
    ExtrapPars = round(ExtrapPars);
end

% Perform extrapolation, if required
if ExtrapolationAllowed
    if T < eph.T(1)
        % Early ephemeris extrapolation
        Extrapolating = true;
        n0 = 1;
        dTextrap  = eph.T(1)-T;
        dTstep    = eph.T(2)-eph.T(1);
    elseif T > eph.T(eph.N)
        % Late ephemeris extrapolation
        Extrapolating = true;
        n0 = eph.N;
        dTextrap  = T-eph.T(eph.N);
        dTstep    = eph.T(eph.N)-eph.T(eph.N-1);
    else
        % No extrapolation required
        Extrapolating = false;
    end
    if Extrapolating
        % Issue extrapolation error or warning
        if dTextrap > ExtrapPars(2)*dTstep
            error(['Attempting extrapolation exceeding maximum allowed ' num2str(ExtrapPars(2)) ' timesteps']);
        elseif dTextrap > ExtrapPars(1)*dTstep
            warning(['Performing extrapolation exceeding ' num2str(ExtrapPars(1)) ' timesteps']);
        end
        % Get eph. time units
        if ~isfield(eph,'SecPerTimeUnit') || isempty(eph.SecPerTimeUnit)
            SecPerTimeUnit = 86400; % Assume time units of days
        else
            SecPerTimeUnit = eph.SecPerTimeUnit;
        end
        % Assume rectilinear, force-free motion for the extrapolation
        r0 = eph.X(1:3,n0)'; % Initial position        
        vi = eph.X(4:6,n0)'; % Constant velocity for force-free motion
        ri = r0 + vi*(T-eph.T(n0))*SecPerTimeUnit;
        % Remaining output parameters
        na = n0; nb = n0; ra = ri; rb = ri; va = vi; vb = vi;
        return;
    end
end

% Find best five bracketing eph points for interpolation
[~,na] = min(abs(T-eph.T));
n1 = na-2;
n2 = na+2;
if (n1 < 1)
    n1 = 1; n2 = 5;
elseif (n2 > eph.N)
    n1 = eph.N-4; n2 = eph.N;
end

% Interpolate state (but not covariance)
i = n1:n2;
[ra,va] = StateCovInterp(T, eph.T(i), eph.X(1:3,i)', eph.X(4:6,i)');

% If appropriate, interpolate states with the alternate five bracketing
% points, and create weighted average states
if NoWeighting || (T == eph.T(na))
    nb = na;
    rb = ra; vb = va;
    ri = ra; vi = va;
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
        % The "alternate" bracketing points are same as the original five
        nb = na;
        rb = ra; vb = va;
        ri = ra; vi = va;
    else
        % Alternate bracketing points are different than the originals,
        % so perform the interpolation
        i = k1:k2;
        [rb,vb] = StateCovInterp(T, eph.T(i), eph.X(1:3,i)', eph.X(4:6,i)');
        % Weight the two alternate sets of states by time differences
        w = (T-eph.T(na)) / (eph.T(nb)-eph.T(na));
        if (w <= 0) || (w >= 1)
            error('Invalid weight');
        end
        omw = 1-w;
        ri = omw*ra+w*rb; vi = omw*va+w*vb;
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
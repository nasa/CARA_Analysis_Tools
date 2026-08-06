function [Xt,out] = convert_ephoffset_to_cartesian(t,EOt,MeanXt,eph,verbose)
% convert_ephoffset_to_cartesian - Convert ephemeris-offset (EO) state to
%                                  cartesian state.
%
% Syntax: [Xt, out] = convert_ephoffset_to_cartesian(t, EOt, MeanXt, eph);
%         [Xt, out] = convert_ephoffset_to_cartesian(t, EOt, MeanXt, eph, verbose);
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
%    t          -   Time (eph.T units, e.g., days)
%
%    EOt        -   Ephemeris offset state                  [6x1]
%                   (eph.X units, e.g., km and km/s)
%
%    MeanXt     -   Mean trajectory cartesian state         [6x1]
%                   (same units as EOt)
%
%    eph            -   Ephemeris structure with fields:
%                         .T - Ephemeris times              [1xN]
%                         .X - State vectors                [6xN]
%                         .P - Covariance matrices          [6x6xN]
%                         .N - Number of points
%
%    verbose    -   Verbosity level: 0, 1, 2, or 3
%                   (optional, default = 0)
%
% =========================================================================
%
% Output:
%
%   Xt          -   Cartesian state transformed from EOt   [6x1]
%
%   out         -   Auxiliary output
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Initializations
Nargin = nargin;
if Nargin < 5; verbose = []; end
if isempty(verbose); verbose = 0; end % Levels 0, 1, 2, or 3

% Persistent parameters
persistent ExtrapPars

% Define the ephemeris extrapolation parameters
if isempty(ExtrapPars)
    % Issue a warning for extrapolations beyond first number of eph.
    % points, but issue no errors because eph. offset coordinates can
    % sometimes require state extrapolations
    ExtrapPars = [Inf Inf];
end

% Get eph. time units
if ~isfield(eph,'SecPerTimeUnit') || isempty(eph.SecPerTimeUnit)
    SecPerTimeUnit = 86400; % Assume time units of days
else
    SecPerTimeUnit = eph.SecPerTimeUnit;
end

% Squared mean velocity at time t
vbart = MeanXt(4:6);
vbartmag2 = vbart' * vbart;
vbartmag = sqrt(vbartmag2);

% Calculate the ephemeris time, tau, from the first EO element
tau = t + EOt(1)/vbartmag/SecPerTimeUnit;

% Interpolate the mean position at tau from ephemeris
[rbartau,vbartau] = interpState(tau,eph,[],ExtrapPars);

% Non-convergence
if isempty(rbartau)
    if verbose > 0
        warning('Offset time out of ephemeris bounds');
    end
    Xt = NaN(6,1); out = [];
    return;
end

% Convert to column vectors
rbartau = rbartau'; vbartau = vbartau';

% Calculate the along-traj. vector (like V of VNB)
Xhat_tau = vbartau/norm(vbartau);
% Calculate first traj.-normal vector (like N of VNB)
Yhat_tau = cross(rbartau,vbartau); Yhat_tau = Yhat_tau/norm(Yhat_tau);
% Calculate second traj.normal vector (like B of VNB)
Zhat_tau = cross(Xhat_tau,Yhat_tau);

% Calculate position vector
rt = rbartau + Yhat_tau*EOt(2) + Zhat_tau*EOt(3);

% Calculate velocity
vt = vbart + Xhat_tau*EOt(4) + Yhat_tau*EOt(5) + Zhat_tau*EOt(6);

% Cartesian state
Xt = [rt; vt];

% Other output
out.tau = tau;
out.Xbar_tau = [rbartau; vbartau];
out.Xhat_tau = Xhat_tau';
out.Yhat_tau = Yhat_tau';
out.Zhat_tau = Zhat_tau';

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
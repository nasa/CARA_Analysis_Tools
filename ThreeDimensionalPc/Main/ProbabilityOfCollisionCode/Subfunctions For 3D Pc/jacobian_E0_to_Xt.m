function [JT,XT] = jacobian_E0_to_Xt(T,E0,fr,mu,errflag)
% =========================================================================
%
% Use epoch equinoctial orbital element state, Ecent0, to calculate the 
% cartesian states, XcentT = (rvec',vvec')', for a set of offset times from
% epoch, T = t-t0.  Here Ecent0 represents the "center-of-expansion"
% equinoctial elements used for the first-order Taylor series expansion of
% the motion:
%              X(t) = Xcent(t) + Jcent(t) * (E0 - Ecent0)
% or
%              E0 = Ecent0 + [Jcent(t)]^-1 * (X(t)-Xcent(t))
% with
%              J(t) = [dX(t)/dE0] at E0 = Ecent0 (i.e., X(t) = Xcent(t))
%
% =========================================================================
%
% INPUT:
%
% T  = Time offsets from initial (s). [NTx1] or [1xNT]
% E0 = Equinoctial elements at initial time t0. [6x1]
%      [n,af,ag,chi,psi,lam]' at t = t0.
%      See Vallado and Alfano (2015) for details.
% fr = Equinoctial element retrograde factor (optional, default = +1)
% mu = Gravitational constant (optional).
% errflag = Error flag (optional, default = 2)
%   0 => No error or warning issued for F nonconvergence (not recommended)
%   1 => Warning issued for F nonconvergence (not recommended)
%   2 => Error   issued for F nonconvergence (recommended)
%
% =========================================================================
%
% OUTPUT:
%
% JT = EpochEquinoctial-to-EphemerisCartesian transformations [6x6xNT]
% XT = Center-of-expansion cartesian state ephemeris (km & km/s) [6xNT]
%
% =========================================================================
%
% REFERENCE:
%
% Vallado and Alfano (2015), AAS 15-537
% Broucke and Cefola (1972), Celestial Mechanics, Vol.5, pp.303-310
%
% =========================================================================

% Defaults and intializations

Nargin = nargin; % Nargout = nargout;

na = 3;
if Nargin < na || isempty(fr)
    % Default to prograde equinoctial elements
    fr = 1;
end

na = na+1;
if Nargin < na || isempty(mu)
    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu  = 3.986004418e5;
end

na = na+1;
if Nargin < na || isempty(errflag)
    errflag = 2;
end

% Number of ephemeris times
NT = numel(T);

% Calculate X(t)
[rT,vT] = convert_equinoctial_to_cartesian( ...
    E0(1),E0(2),E0(3),E0(4),E0(5),E0(6),T,  ...
    fr,mu,[],[],errflag);
XT = [rT; vT];

% Initialize dE(t)/dE(t0) STM
phi = eye(6,6);

% Initialize E(t)
ET = E0;

% Initialize output Jacobian array
JT = NaN(6,6,NT);

% Loop over times and calculate output

for nT=1:NT
    % Calculate the equinoctial mean longitude at this time
    ET(6) = E0(6)+T(nT)*E0(1);
    % Define off-diagonal dE(t)/dE(t0) STM element
    phi(6,1) = T(nT);
    % Calculate dX(t)/dE(t) STM
    J = jacobian_equinoctial_to_cartesian(ET,XT(:,nT),fr,mu);
    % Calculate the dX(t)/dE0 STM
    JT(:,:,nT) = J*phi;
end

return
end
function [a,n,af,ag,chi,psi,lM,F] = ...
    convert_cartesian_to_equinoctial(rvec,vvec,fr,mu,issue_warnings)
% =========================================================================
%
% Convert cartesian state (r,v) to the equinoctial orbital elements 
% (a,n,af,ag,chi,psi,lambdaM,F)
%
% =========================================================================
%
% INPUT:
%
% rvec = Position vector (km) [3x1]
% vvec = Velocity vector (km] [3x1]
% fr   = Equinoctial element retrograde factor (optional, default = +1)
% mu   = Gravitational constant (optional)
%
% =========================================================================
%
% OUTPUT:
%
% (a,n,af,ag,chi,psi,lambdaM,F) = Equinoctial elements and related
% quantities.  See Vallado and Alfano (2015) for details.
%
% =========================================================================
%
% REFERENCE:
%
% Vallado and Alfano (2015), AAS 15-537
%
% =========================================================================

% Defaults and intializations

Nargin = nargin;

if Nargin < 3 || isempty(fr)
    % Default to prograde equinoctial elements
    fr = 1;
end

if Nargin < 4 || isempty(mu)
    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu = 3.986004418e5;
end

if Nargin < 5 || isempty(issue_warnings)
    issue_warnings = true;
end

% Place Vectors in correct format
rvec = reshape(rvec,3,1);
vvec = reshape(vvec,3,1);

% Calculate equinoctial elements using equations for bound orbits.
% See Vallado and Alfano (2015) equations (9) on page 6.

r2 = rvec'*rvec;
r = sqrt(r2);

v2 = vvec'*vvec;
% v = sqrt(v2);

a = mu*r/(2*mu-v2*r);

% Issue warning and return for unbound orbits: a <= 0 implies E >= 0
if (a <= 0)
    if issue_warnings
        warning('Cannot process unbound orbit with nonpositive semimajor axis (i.e., a <= 0)');
    end
    n=[]; af=[]; ag=[]; chi=[]; psi=[]; lM=[]; F=[];
    return;
end

rdv = rvec' * vvec;
rcv = cross(rvec,vvec);

% Calculating remaining equinoctial elements for this bound orbit
    
n = sqrt(mu/a^3);

evec = (1/mu) * ( (v2-mu/r)*rvec - rdv*vvec );

what = rcv/norm(rcv);
we = what(1); 
wq = what(2);
ww = what(3);

cpden = (1+fr*ww);

chi =  we/cpden;
psi = -wq/cpden;

chi2 = chi^2; psi2 = psi^2; C = 1 + chi2 + psi2;

fhat = [ 1-chi2+psi2  ; 2*chi*psi        ; -2*fr*chi ] / C;
ghat = [ 2*fr*chi*psi ; (1+chi2-psi2)*fr ; 2*psi     ] / C;

% [fhat ghat what]

af = fhat' * evec;
ag = ghat' * evec;

af2 = af^2; ag2 = ag^2;

X = fhat' * rvec;
Y = ghat' * rvec;

safg = sqrt(1-ag2-af2);
b = 1 / (1+safg);
Fden = a*safg;

bagaf = b*ag*af;

sinF = ag + ((1-ag2*b)*Y-bagaf*X)/Fden;
cosF = af + ((1-af2*b)*X-bagaf*Y)/Fden;
F = atan2(sinF,cosF);

lM = F + ag*cosF - af*sinF;

return
end
function [a,n,af,ag,chi,psi,lM,F] = ...
    convert_cartesian_to_equinoctial(rvec,vvec,fr,mu,issue_warnings)
% convert_cartesian_to_equinoctial - Convert cartesian state (r,v) to the 
% equinoctial orbital elements (a,n,af,ag,chi,psi,lambdaM,F)
%
% Syntax: [a,n,af,ag,chi,psi,lM,F] = ...
%         convert_cartesian_to_equinoctial(rvec,vvec,fr,mu,issue_warnings);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   rvec            -   Cartesian position vector (km)                [3x1]
%   vvec            -   Cartesian velocity vector (km)                [3x1]
%   fr              -   Equinoctial element retrograde factor (optional, 
%                       default = +1)
%   mu              -   Gravitational constant (optional, default = 
%                       3.986004418e5)
%   issue_warnings  -   Display warning messages (optional, default = true)
%
% =========================================================================
%
% Output:
%
%   (a,n,af,ag,chi,psi,lambdaM,F) = Equinoctial elements and related
%   quantities. See Vallado and Alfano (2015) for details.  
%   
% =========================================================================
%
% References:
%
%   Vallado and Alfano (2015), AAS 15-537
%
% =========================================================================
%
% Initial version: Dec 2019;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

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

if abs(fr) ~= 1
    if issue_warnings
        warning('convert_cartesian_to_equinoctial:InvalidRetrogradeFactor', 'fr must be either +1, -1, or not passed in (defaullt = +1)');
    end
    a=[]; n=[]; af=[]; ag=[]; chi=[]; psi=[]; lM=[]; F=[];
    return;
end

% Place vectors in column vector format
if ~isequal(size(rvec),[3 1])
    rvec = reshape(rvec,3,1);
end

if ~isequal(size(vvec),[3 1])
    vvec = reshape(vvec,3,1);
end

% Calculate equinoctial elements using equations for bound orbits.
% See Vallado and Alfano (2015) equations (9) on page 6.

r = sqrt(rvec(1)^2 + rvec(2)^2 + rvec(3)^2);

v2 = vvec(1)^2 + vvec(2)^2 + vvec(3)^2;

a = mu*r/(2*mu-v2*r);

% Issue warning and return for unbound orbits: a <= 0 implies E >= 0
if (a <= 1E-5 || a == Inf)
    if issue_warnings
        warning('convert_cartesian_to_equinoctial:UnboundOrbit', 'Cannot process unbound orbit with infinite or nonpositive semimajor axis (a)');
    end
    n=[]; af=[]; ag=[]; chi=[]; psi=[]; lM=[]; F=[];
    return;
end

rdv = rvec(1) * vvec(1) + rvec(2) * vvec(2) + rvec(3) * vvec(3);

% rxv cross product
% rcv = cross(rvec,vvec); % slower way
rcv = [rvec(2)*vvec(3)-rvec(3)*vvec(2); ...
       rvec(3)*vvec(1)-rvec(1)*vvec(3); ...
       rvec(1)*vvec(2)-rvec(2)*vvec(1)];

% Calculating remaining equinoctial elements for this bound orbit
    
n = sqrt(mu/a^3);

evec = (1/mu) * ( (v2-mu/r)*rvec - rdv*vvec );

% Issue warning and return for unbound orbits with ecc^2 >= 1
if evec'*evec >= 1
    if issue_warnings
        warning('convert_cartesian_to_equinoctial:UnboundOrbit', 'Cannot process unbound orbit with ecc^2 = evec.evec >= 1');
    end
    n=[]; af=[]; ag=[]; chi=[]; psi=[]; lM=[]; F=[];
    return;
end

% Calculate chi and psi elements
what = rcv/sqrt(rcv(1)^2 + rcv(2)^2 + rcv(3)^2);

if what(3) + fr <= 1E-10
    if issue_warnings
        warning('convert_cartesian_to_equinoctial:EquatorialRetrogradeOrbit', 'Cannot process equatorial retrograde orbits (i = 180) without flag');
    end
    n=[]; af=[]; ag=[]; chi=[]; psi=[]; lM=[]; F=[];
    return;
end

cpden = (1+fr*what(3));
chi =  what(1)/cpden;
psi = -what(2)/cpden;

chi2 = chi^2; psi2 = psi^2; C = 1 + chi2 + psi2;

% Calculte f and g unit vectors
fhat = [ 1-chi2+psi2  ; 2*chi*psi        ; -2*fr*chi ] / C;
ghat = [ 2*fr*chi*psi ; (1+chi2-psi2)*fr ; 2*psi     ] / C;

% [fhat ghat what]

af = fhat(1) * evec(1) + fhat(2) * evec(2) + fhat(3) * evec(3);
ag = ghat(1) * evec(1) + ghat(2) * evec(2) + ghat(3) * evec(3);

af2 = af^2; ag2 = ag^2; ec2 = ag2+af2;

% Issue warning and return for unbound orbits with ecc^2 >= 1
if ec2 >= 1
    if issue_warnings
        warning('convert_cartesian_to_equinoctial:UnboundOrbit', 'Cannot process unbound orbit with ecc^2 = af^2 + ag^2 >= 1');
    end
    n=[]; af=[]; ag=[]; chi=[]; psi=[]; lM=[]; F=[];
    return;
end

X = fhat(1) * rvec(1) + fhat(2) * rvec(2) + fhat(3) * rvec(3);
Y = ghat(1) * rvec(1) + ghat(2) * rvec(2) + ghat(3) * rvec(3);

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

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 12-13-2019 | Initial development
% D. Hall        | 07-29-2020 | Added issue_warnings parameter
% D. Hall        | 04-12-2022 | Changed some calculations for speed
% E. White       | 06-30-2023 | Added compliant documentation, added input
%                               validation, improved check for unbound
%                               orbits, improved speed in various
%                               locations, added warning for retrograde
%                               orbits without flag
% D. Hall        | 02-27-2024 | Added processing to catch unbound orbits
%                               with ecc >= 1, and prevent errors for
%                               orbits with e^2 values within round-off
%                               error of 1.

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

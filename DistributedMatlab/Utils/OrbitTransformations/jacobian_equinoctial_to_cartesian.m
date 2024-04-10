function J = jacobian_equinoctial_to_cartesian(E,X,fr,mu)
% =========================================================================
%
% Calculate the Jacobian, J = dX/dE, between the equinoctial elements
%
%   E = (n,af,ag,chi,psi,lM)'   [ Order used by Vallado & Alfano (2015) ]
%
% and cartesian state
%
%   X = (r',v')'                [ usually km & km/s, unless mu is input in
%                                 different units ]
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% INPUT:
%
% E = (n,af,ag,chi,psi,lM)' = Equinoctial state vector (km) [6x1]
%                             (The mean motion n has units of radians/s)
% X = (r',v')' = Cartesian state vector (km, km/s) [6x1]
% fr   = Equinoctial element retrograde factor (optional, default = +1)
% mu   = Gravitational constant (optional)
%
% =========================================================================
%
% OUTPUT:
%
% J = dX/dE [6x6] = Jacobian as summarized by Vallado & Alfano (2015), 
% appendix "Equinoctial to Cartesian" section, p.23
%
% =========================================================================
%
% REFERENCE:
%
% Vallado and Alfano, "Updated Analytical Partials for Covariance
%     Transformations and Optimization, 2015, AAS 15-537
%
% =========================================================================
%
% Initial version: Dec 2019;  Latest update: Mar 2024
%
% ----------------- BEGIN CODE -----------------

% Defaults and intializations

Nargin = nargin;

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

% Place Vectors in correct format
X = reshape(X,6,1);

% Extract the elements

n = E(1); af = E(2); ag = E(3); chi = E(4); psi = E(5); % lM = E(6);

rvec = X(1:3); vvec = X(4:6);

r2 = rvec'*rvec; r = sqrt(r2); r3 = r2*r;

% Aux quantities

n2 = n^2; a3 = mu/n2; a = nthroot(a3,3); A = n*a^2;

ag2 = ag^2; af2 = af^2; B = sqrt(1-ag2-af2); % b = 1/(1+B);

chi2 = chi^2; psi2 = psi^2; C = 1 + chi2 + psi2;

% fhat and ghat

fhat = [ 1-chi2+psi2  ; 2*chi*psi        ; -2*fr*chi        ] / C;
ghat = [ 2*fr*chi*psi ; (1+chi2-psi2)*fr ; 2*psi            ] / C;
what = [ 2*chi        ; -2*psi           ; (1-chi2-psi2)*fr ] / C;

% X, Y, Xdot, and Ydot

X = fhat' * rvec;
Y = ghat' * rvec;

Xd = fhat' * vvec;
Yd = ghat' * vvec;

% Partials of X, Y, Xdot, Ydot

AB = A*B; Bp1 = B+1; nBp1 = n*Bp1; Aor3 = A/r3;

dXdaf =  ag*Xd/nBp1 + a*(Y*Xd/AB - 1);
dYdaf =  ag*Yd/nBp1 - a*(X*Xd/AB);
dXdag = -af*Xd/nBp1 + a*(Y*Yd/AB);
dYdag = -af*Yd/nBp1 - a*(X*Yd/AB + 1);

dXddaf =  a*Xd*Yd/AB - Aor3*(a*ag*X/Bp1+X*Y/B);
dYddaf = -a*Xd*Xd/AB - Aor3*(a*ag*Y/Bp1-X*X/B);
dXddag =  a*Yd*Yd/AB + Aor3*(a*af*X/Bp1-Y*Y/B);
dYddag = -a*Xd*Yd/AB + Aor3*(a*af*Y/Bp1+X*Y/B);

% Define the Jacobian, J = dX/dE

J = NaN(6,6);

% drvec/dn & dvvec/dn
cv = 1/(3*n);
cr = -2*cv;
J(1:3,1) = cr*rvec;
J(4:6,1) = cv*vvec;

% drvec/daf & dvvec/daf
J(1:3,2) = dXdaf *fhat + dYdaf *ghat;
J(4:6,2) = dXddaf*fhat + dYddaf*ghat;

% drvec/dag & dvvec/dag
J(1:3,3) = dXdag *fhat + dYdag *ghat;
J(4:6,3) = dXddag*fhat + dYddag*ghat;

% drvec/dchi & dvvec/dchi
cc = 2/C;
J(1:3,4) = cc * ( fr*psi*(Y *fhat-X *ghat) - X *what );
J(4:6,4) = cc * ( fr*psi*(Yd*fhat-Xd*ghat) - Xd*what );

% drvec/dpsi & dvvec/dpsi
J(1:3,5) = cc * ( fr*chi*(X *ghat-Y *fhat) + Y *what );
J(4:6,5) = cc * ( fr*chi*(Xd*ghat-Yd*fhat) + Yd*what );

% drvec/dlambdaM & dvvec/dlambdaM
J(1:3,6) = vvec/n;
J(4:6,6) = (-n*a3/r3)*rvec ;

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
% L. Baars       | 03-05-2024 | Modified fractional exponent to use more
%                               stable nthroot() function.
% L. Baars       | 03-27-2024 | Modified header/footer to standard
%                               documentation and made some minor comment
%                               changes.

% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
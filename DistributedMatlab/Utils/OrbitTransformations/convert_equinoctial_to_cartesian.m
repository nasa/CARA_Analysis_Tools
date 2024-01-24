function [rvec,vvec,X,Y,Xdot,Ydot,F,cF,sF] = ...
    convert_equinoctial_to_cartesian(n,af,ag,chi,psi,lam0,T, ...
                                     fr,mu,Ftol,maxiter,errflag)
% convert_equinoctial_to_cartesian - Convert equinoctial orbital elements 
% (n,af,ag,chi,psi,lam0) to cartesian states (rvec,vvec) for a set of 
% offset times from epoch (T).
%
% Syntax: [rvec,vvec,X,Y,Xdot,Ydot,F,cF,sF] = ...
%         convert_equinoctial_to_cartesian(n,af,ag,chi,psi,lam0,T,fr, ...
%         mu,Ftol,maxiter,errflag);
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
%   (n,af,ag,chi,psi,lam0)  -   Epoch equinoctial elements (see       [1x1]
%                               Vallado and Alfano (2015) for details) 
%   T                       -   Time offsets from epoch (s)  [nx1] or [1xn]
%   fr                      -   Equinoctial element retrograde factor 
%                               (optional, default = +1)
%   mu                      -   Gravitational constant (optional, default =
%                               3.986004418e5)
%   Ftol                    -   Tolerance for F angle convergence 
%                               (optional, default = 100*eps(2*pi))
%   maxiter                 -   Max iterations for F convergence (optional, 
%                               default = 100)
%   errflag                 -   Error flag (optional, default = 2)
%                                   0 => No error or warning issued for F 
%                                        nonconvergence (not recommended)
%                                   1 => Warning issued for F 
%                                        nonconvergence (not recommended)
%                                   2 => Error issued for F nonconvergence 
%                                        (recommended)
%
% =========================================================================
%
% Output:
%
%   rvec                    -   Cartesian position vector             [3xn]  
%   vvec                    -   Cartesian velocity vector             [3xn]
%   (X,Y,Xdot,Ydot,F,cF,sF) -   Quantities from the calculation       [1xn]
%                               (See Vallado and Alfano (2015) for details)
% =========================================================================
% 
% References (optional):
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

na = 8;
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
if Nargin < na || isempty(Ftol)
    Ftol = 100*eps(2*pi);
end

na = na+1;
if Nargin < na || isempty(maxiter)
    maxiter = 100;
end

na = na+1;
if Nargin < na || isempty(errflag)
    errflag = 2;
end

% Semimajor axis and related quantities

n2 = n.^2;
a3 = mu./n2;
a = a3.^(1/3);
na = n.*a;

% Mean longitudes

lam = lam0 + n.*T;

% (ag,af) quantities

ag2 = ag.^2; af2 = af.^2;

B = (1-ag2-af2).^0.5;
b = 1./(1+B);

omag2b = 1-ag2.*b;
omaf2b = 1-af2.*b;
afagb = af.*ag.*b;

% Unit vectors (f,g)

chi2 = chi.^2; psi2 = psi.^2; C = 1+chi2+psi2;

fhat = [ 1-chi2+psi2  ; 2.*chi.*psi        ; -2*fr.*chi ] ./ C;
ghat = [ 2*fr.*chi.*psi ; (1+chi2-psi2)*fr ; 2.*psi     ] ./ C;

% Calculate states and output quantities for each ephemeris point

NT = numel(T);

rvec = NaN(3,NT);
vvec = NaN(3,NT);

X = NaN(1,NT); Y = X; Xdot = X; Ydot = X;
F = X; cF = X; sF = X;

for nt=1:NT
    
    % Eccentric longitude
    
    [F(nt),~,cF(nt),sF(nt)] = ...
        equinoctial_kepeq(lam(nt),af,ag,Ftol,maxiter,errflag);
    
    X(nt) = a * (omag2b.*cF(nt) + afagb.*sF(nt) - af);
    Y(nt) = a * (omaf2b.*sF(nt) + afagb.*cF(nt) - ag); 
    
    rvec(1:3,nt) = X(nt)*fhat + Y(nt)*ghat;

    na2or = na/(1-af*cF(nt)-ag*sF(nt));

    Xdot(nt) = na2or * (afagb*cF(nt) - omag2b*sF(nt));
    Ydot(nt) = na2or * (omaf2b*cF(nt) - afagb*sF(nt));
    
    vvec(1:3,nt) = Xdot(nt)*fhat + Ydot(nt)*ghat;

end

return
end

function [F,converged,cF,sF] = equinoctial_kepeq(lam,af,ag,Ftol,maxiter,errflag)

% Solve Kepler's equation in equinoctial elements.

% Ensure mean longitude in the range 0 <= lam < 2*pi
lam = mod(lam,2*pi);

% Initialize convergence flag
converged = false;

% Calculate eccentricity squared
ecc2 = af^2+ag^2;

if ecc2 >= 1
    % Handle warning or error
    errtext = 'equinoctial_kepeq cannot process eccentricities >= 1' ;
    if errflag == 1
        warning(errtext);
        F = []; cF = []; sF = [];
        return;
    elseif errflag == 2
        error(errtext);
    elseif errflag ~= 0
        error('Invalid errflag value');
    end
end

% if ecc2 < 0.81 % Corresponds to ecc < 0.9 (used for testing)
    
    % Use iterative algorithm of Valado and Alfano (2015) equation 10 
    % to find eccentric longitude F

    F = lam; cF = cos(F); sF = sin(F);
    
    iter = 0;
    still_looping = true;

    while still_looping

        % Calculate the increment for F
        % top  = F + ag*cF - af*sF - lam;
        % bot  = 1 - ag*sF - af*cF;
        % Fdel = top/bot;
        Fdel = (F + ag*cF - af*sF - lam) / (1 - ag*sF - af*cF);

        % Update F
        F = F - Fdel; cF = cos(F); sF = sin(F);

        % Check for F convergence
        absFdel = abs(Fdel);
        if (absFdel < Ftol)
            % Nominal convergence
            converged = true; still_looping = false;
        elseif (iter >= maxiter) || (absFdel >= pi)
            % Failure to converge from to many iterations or for too large
            % an increment
            still_looping = false;
        else
            % Increment iteration counter
            iter = iter+1;
        end
        
        % if iter == 0
        %     maxabsFdel = absFdel;
        % else
        %     maxabsFdel = max(absFdel,maxabsFdel);
        % end

    end

% end

% If solution remains unconverged, then use Matlab's fminbnd function, 
% which is more stable for high eccentricites than algorithm of 
% Valado and Alfano (2015).

if ~converged
    
    % Anonymous function for Kepler's equation in equinoctial elements    
    KEP = @(FF) abs(FF + ag*cos(FF) - af*sin(FF) - lam);

    % Solve Kepler's equation using Matlab's fminbnd function
    options = optimset('fminbnd');
    options = optimset(options,'TolX',Ftol,'MaxIter',max(maxiter,500));
    [F0,~,exitflag] = fminbnd(KEP,0,2*pi,options);
    
    % Check for convergence
    if exitflag == 1
        converged = true;
        F = F0; cF = cos(F); sF = sin(F);
    end
    
    % If solution still remains unconverged, then issue warning or error
    if ~converged
        % Handle warning or error
        errtext = 'equinoctial_kepeq failed to converge' ;
        if errflag == 1
            warning(errtext);
            F = []; cF = []; sF = [];
        elseif errflag == 2
            error(errtext);
        elseif errflag ~= 0
            error('Invalid errflag value');
        end
    end

end

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
% T. Lechtenberg | 03-16-2020 | Modified error handling
% D. Hall        | 07-30-2020 | Added error checking and expanded upon
%                               Kepler's equation
% D. Hall        | 04-12-2022 | Formatted and reorganized code
% E. White       | 06-30-2023 | Added compliant documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

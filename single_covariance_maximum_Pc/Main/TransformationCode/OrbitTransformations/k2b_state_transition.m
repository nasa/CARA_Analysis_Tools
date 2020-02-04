function [params,converged,x,Phi] = k2b_state_transition(x0,T,mu,params)
% =========================================================================
%
% Given an initial state, x0, calculate an ephemeris for Keplerian
% 2-body motion, which includes the states x and, optionally, the error 
% or state-deviation state-transition matrices, Phi.
%
% Algorithm taken from Shepperd (1985) - see reference below.
%
% All units MKS unless otherwise specified.
%
% =========================================================================
%
% INPUT:
%
%   x0      = State vector at time = t0                     [6 x 1].
%   T       = Time-of-flight (TOF) from time = t0           [1 x NT].
%   mu      = Grav. constant times central-body mass        [1 x 1].
%   params  = Structure holding parameters and (optionally)
%             parts of a previous calculation. Using this
%             can speed repetitive calculations that use
%             the same x0 vector.
%
% =========================================================================
%
% OUTPUT:
%
%   params  = Updated params structure.
%   x       = State vectors at times = t+t0                 [6 x NT].
%   Phi     = State transition matrices at times = t+t0     [6 x 6 x NT].
%
% =========================================================================
%
% REFERENCES:
%
%    S. W. Shepperd, "Universal Keplerian State Transition Matrix"
%    Celestial Mechanics, 35, 129-144, 1985.  [Hereafter S85]
%
%    W. H. Goodyear, "Completely General Closed-Form Solution for
%    Coordinates and Partial Derivatives of the Two-Body Problem"
%    Astron. J., 70, 189-192, 1965.           [Hereafter G65]
%
% =========================================================================
%
%
% =========================================================================

% Check to see if the params structure needs to be updated. Specifically,
% if the input params structure was originally built using the same state,
% x0, then many of the quantities saved in params do not have to be
% recalculated, thereby saving CPU time.

if nargin < 4; params = []; end

if isfield(params,'x0') && isequal(params.x0,x0)
    update_params = false;
else
    update_params = true;
end

% Define and/or update params, if required

if update_params
    
    % disp('Updating k2b_state_transition parameters');
    
    % Save the state vector at t0
    
    params.x0 = x0;
    
    % Defaults and constants
    
    params.kepler_max_iterations = 25;
    params.time_rel_tolerance = 1e-12;
    
    params.cont_frac_max_iterations = 50;
    params.cont_frac_tolerance = 1e-9;
    
    params.twopi = 2*pi;
    params.c16o15 = 16/15;

    % Extract pos and vel vectors at t0

    params.r0vec = x0(1:3);
    params.r02 = params.r0vec' * params.r0vec;
    params.r0 = sqrt(params.r02);
    
    params.v0vec = x0(4:6);
    params.v02 = params.v0vec' * params.v0vec;
    params.v0 = sqrt(params.v02);

    params.r0v0T = [params.r0vec params.v0vec]';
    
    % S85 eqs (A.2) & (A.3)
    
    params.nu0 = params.r0vec' * params.v0vec;
    
    params.beta = 2*mu/params.r0 - params.v02;
    
    % S85 eq (A.6) for period, and mean motion
    
    params.P = params.twopi * mu *(params.beta)^(-1.5);
    params.a0 = mu / params.r0^3;
    
    params.time_tolerance = params.time_rel_tolerance*params.P;
    
    % S85 eq (A.7), 2nd + 3rd terms that do not depend on T
    
    params.n2 = 0.5 - 2*params.nu0/(params.beta*params.P);
    
    % Max/min u
    
    if (params.beta ~= 0)
       params.umax = 1 / sqrt(abs(params.beta));
    else
       params.umax = 1.0e24;
    end
    params.umin = -params.umax;
    
    % Misc
    
    params.I3x3 = eye(3,3);
    
    params.twopi_beta_m2p5 = params.twopi * params.beta^(-2.5);
    
end

% Number of ephemeris time points. These are times-of-flight relative to
% t0, the time when the state is x0.

NT = numel(T);

% Initialize output states.

x = zeros(6,NT);

% Initialize convergence flags

converged = false(size(T));

% Check to see if state transition matrices are desired.  If so initialize
% the required arrays.

calc_Phi = (nargout > 3);
if calc_Phi
    Phi = zeros(6,6,NT);
    M = zeros(3,3);
end

% Calculate S85 DeltaU (A.8) for all of the times

if (params.beta <= 0)
    DeltaU = zeros(size(T));
else
    n = floor(T./params.P + params.n2);
    DeltaU = params.twopi_beta_m2p5 * n;
end

% Process all of the times

for nt=1:NT
    
    % Current time
    
    Tnow = T(nt);
    
    % Initialization S85 (A.4)
    
    u = 0;
    
    % Initialize quantities derived from u for the first iteration.
    % This saves having to calculate k2b_cont_frac(q=0) = 1 and wasting
    % CPU time.
    
    q = 0;            % q = params.beta * u^2;
    U0_halfw = 1;     % U0_halfw = 1 - 2 * q;
    U1_halfw = 0;     % U1_halfw = 2 * u * (1 - q);
    G = 1;            % G = k2b_cont_frac(q=0);    
    
    % Initialize u limits from saved parameters.
    
    umax = params.umax;
    umin = params.umin;
    
    % Initialize flag for continued fraction iteration convergence
    
    cont_frac_converged = true;

    % Kepler iteration loop, S85 (A.9)-(A.19)
    
    for niter=1:params.kepler_max_iterations
        
        % S85 (A.12)
        
        U = params.c16o15 * U1_halfw^5 * G + DeltaU(nt);
        
        % S85 (A13)-(A15)
        
        U0 = 2 * U0_halfw^2 - 1;
        U1 = 2 * U0_halfw * U1_halfw;
        U2 = 2 * U1_halfw^2;
        U3 = params.beta*U + U1*U2/3;

        % S85 (A17)-(A.18)

        r = params.r0*U0 + params.nu0*U1 + mu*U2;
        t = params.r0*U1 + params.nu0*U2 + mu*U3;
        
        % Check for time convergence

        dt = Tnow-t;
        
        % disp(['niter = ' num2str(niter) ' T = ' num2str(Tnow) ' t = ' num2str(t) ' dt = ' num2str(dt)]);
        
        if abs(dt) < params.time_tolerance
            
            % If time has converged, then set the Kepler iteration 
            % convergence flag equal to the continued fraction convergence
            % flag
            
            converged(nt) = cont_frac_converged;
            break;
            
        else
            
            % S85 (A.19)

            du = dt/(4*(1-q)*r);
            
            % % Update u, using simple addition
            % 
            % u = u+du;

            % Update u, but prevent violating the min/max boundaries
            
            if (du < 0)
                umax = u;
                u = u + du;
                if (u < umin)
                    u = 0.5 * (umin + umax);
                end   
            else
                umin = u;
                u = u + du;
                if (u > umax)
                    u = 0.5 * (umin + umax);
                end   
            end
            
            % S85 (A.9)

            q = params.beta * u^2;
            q = q / (1 + q);
            
            % S85 (A.10)-(A.11)

            U0_halfw = 1 - 2 * q;
            U1_halfw = 2 * u * (1 - q);

            % S85 continued fraction 

            [G,cfconv] = k2b_cont_frac(q, ...
                    params.cont_frac_tolerance, ...
                    params.cont_frac_max_iterations);
                
            % See if continue fraction iterations converged. If not, then
            % set the overall convergence flag to false.
                
            if ~cfconv
                cont_frac_converged = false;
            end

        end

    end

    % Calculate the Kepler 2-body state solution using the (f,g) functions,
    % S85 (A.34)-(A.39)

    f = 1 - (mu/params.r0)*U2;
    g = params.r0*U1 + params.nu0*U2;
    F = -mu*U1/(r*params.r0);
    G = 1 - (mu/r)*U2;
    
    x(1:3,nt) = f*params.r0vec + g*params.v0vec;
    x(4:6,nt) = F*params.r0vec + G*params.v0vec;
    
    % Calculate the state transition matrix if required
    
    if calc_Phi
        
        % S85 (A.40)-(A.46)
        
        W = g*U2 + 3*mu*U;
        
        a1 = mu / r^3;
        
        fm1 = f-1;
        Gm1 = G-1;

        % M(1,1) = F * (U0 / (params.r0 * r) + 1 / params.r02 + 1 / (r * r));
        M(1,1) = F * ( ( U0/params.r0 + 1/r ) / r + 1/params.r02 );
        M(1,2) = (F * U1 + (Gm1 / r)) / r;
        M(1,3) = Gm1 * U1 / r;
        M(2,1) = -(F * U1 + (fm1 / params.r0)) / params.r0;
        M(2,2) = -F * U2;
        M(2,3) = -Gm1 * U2;
        M(3,1) = fm1 * U1 / params.r0;
        M(3,2) = fm1 * U2;
        M(3,3) = g * U2;

        a1W = a1*W;
        M(1,1) = M(1,1) - params.a0 * a1W;
        M(1,3) = M(1,3) - a1W;
        M(3,1) = M(3,1) - params.a0 * W;
        M(3,3) = M(3,3) - W;
        
        rv = [x(1:3,nt) x(4:6,nt)];
        
        % Note, for these operations "rv" is a 3x2 matrix, and
        % "params.r0v0T" is a 2x3 matrix, see S85 (13)-(19) and 
        % associated discussion. Also see G65 for background.

        p11 = f * params.I3x3 + rv * M(2:3,1:2) * params.r0v0T;
        p12 = g * params.I3x3 + rv * M(2:3,2:3) * params.r0v0T;
        p21 = F * params.I3x3 - rv * M(1:2,1:2) * params.r0v0T;
        p22 = G * params.I3x3 - rv * M(1:2,2:3) * params.r0v0T;
        
        % Construct 6x6 state transition matrix for this ephemeris point

        Phi(:,:,nt) = [p11 p12; p21 p22];
        
    end
    
end

if any(~converged)
    error('Unconverged K2B state transition case(s).');
end
    
return;
end
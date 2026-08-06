function [conv,rpk,v1pk,v2pk,aux] = PeakOverlapPosEph(t,xb1,Pb1,xb2,Pb2,eph1,eph2,HBR,params)
% PeakOverlapPosEph - Find the inertial-frame position, rpk, of the
%                     point of peak overlap of the primary and secondary
%                     distributions for time = t, along with the mean
%                     primary and secondary velocities at that point.
%
% Syntax: [conv, rpk, v1pk, v2pk, aux] = PeakOverlapPosEph(t, xb1, Pb1, xb2, Pb2, eph1, eph2, HBR);
%         [conv, rpk, v1pk, v2pk, aux] = PeakOverlapPosEph(t, xb1, Pb1, xb2, Pb2, eph1, eph2, HBR, params);
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
%    t          -   Ephemeris time
%
%    xb1        -   Primary cartesian state                 [6x1]
%                   (km and km/s)
%
%    Pb1        -   Primary covariance matrix               [6x6]
%
%    xb2        -   Secondary cartesian state               [6x1]
%                   (km and km/s)
%
%    Pb2        -   Secondary covariance matrix             [6x6]
%
%    eph1       -   Primary ephemeris structure (empty for
%                   equinoctial element mode)
%
%    eph2       -   Secondary ephemeris structure (empty for
%                   equinoctial element mode)
%
%    HBR        -   Combined hard body radius (km)
%
%    params     -   Parameter structure with optional fields:
%                     .Jb1, .Eb1, .Qb1 - Primary equinoctial
%                       Jacobian, elements, and covariance
%                       (required for equinoctial mode)
%                     .Jb2, .Eb2, .Qb2 - Secondary equinoctial
%                       Jacobian, elements, and covariance
%                       (required for equinoctial mode)
%                     .MD2tol   - Convergence tolerance
%                     .maxiter  - Maximum iterations
%                     .Fclip    - Eigenvalue clipping factor
%                     .GM       - Gravitational constant
%                     .verbose  - Verbosity flag
%                     .UseFasterJacobianAlg - Jacobian algorithm flag
%
% =========================================================================
%
% Output:
%
%   conv        -   Convergence flag (true/false/0.5 for oscillating)
%
%   rpk         -   Peak overlap position                   [3x1]
%
%   v1pk        -   Primary velocity at peak position       [3x1]
%
%   v2pk        -   Secondary velocity at peak position     [3x1]
%
%   aux         -   Auxiliary output structure with fields:
%                     .converged  - Convergence flag
%                     .iteration  - Number of iterations
%                     .failure    - Failure code
%                     .xs1, .xs2  - Expansion-center states
%                     .Js1, .Js2  - Expansion-center Jacobians
%                     .Es1, .Es2  - Expansion-center elements
%                     .Ps1, .Ps2  - Pos/vel covariances
%                     .Sigp       - Peak overlap covariance
%                     .Sigpinv    - Inverse of Sigp
%                     .SigpRem    - Remediated Sigp
%                     .SigpReminv - Inverse of SigpRem
%                     .xu1, .xu2  - Offset states
%                     .dE1, .dE2  - Element differences
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Initializations and defaults
Nargin = nargin; Nargout = nargout;

% Set up default parameters
if Nargin < 9; params = []; end

% Determine if equinoctial (EQ) element or ephemeris offset (EO)
% coordinate element state representation is to be used
UseEphOffsetElements = ~isempty(eph1);
if UseEphOffsetElements && isempty(eph2)
    error('Both primary and secondary ephemeris tables required');
end
UseEquinoctialElements = ~UseEphOffsetElements;

if UseEquinoctialElements
    % If using EQ elements, then exact other info from parameters
    Jb1 = params.Jb1;
    Eb1 = params.Eb1;
    Qb1 = params.Qb1;
    Jb2 = params.Jb2;
    Eb2 = params.Eb2;
    Qb2 = params.Qb2;   
end

% Tolerances for convergence using relative Maha. distance squared
if ~isfield(params,'MD2tol') || isempty(params.MD2tol)
    params.MD2tol = [1e-6 1e-3];
end
if numel(params.MD2tol) == 1
    params.MD2tol = [params.MD2tol min(sqrt(params.MD2tol),3e-2)];
end

% Maximum number of iterations to perform
if ~isfield(params,'maxiter') || isempty(params.maxiter)
    params.maxiter = 100;
end
avgiter = min(35,round(params.maxiter*0.35));
acciter = min(25,round(params.maxiter*0.25));
osciter = min(15,round(params.maxiter*0.15));

% Matrix to remediate SIGMAp maxtrices
SigpRem0 = diag(repmat(HBR,[1 3]).^2);

% Eigenvalue clipping factor
if ~isfield(params,'Fclip') || isempty(params.Fclip)
    params.Fclip = 1e-4;
end
Lclip = (HBR*params.Fclip)^2;

if ~isfield(params,'GM'); params.GM = []; end
if isempty(params.GM)
    % Earth gravitational constant mu = GM (EGM-96) [km^3/s^2]
    params.GM = 3.986004418e5;
end
GM = params.GM;

% Verbosity
if ~isfield(params,'verbose') || isempty(params.verbose)
    params.verbose = false;
end
verbose = params.verbose;

% Other initializations
twopi = 2*pi;
I3x3 = eye(3,3);
Z3x3 = zeros(3,3);

% Jacobian algorithm
if ~isfield(params,'UseFasterJacobianAlg') || ...
   isempty(params.UseFasterJacobianAlg)
    params.UseFasterJacobianAlg = true;
end
UseFasterJacobianAlg = params.UseFasterJacobianAlg;

if UseEquinoctialElements

    % Sin and Cos values for initial mean longitudes
    sinLb1 = sin(Eb1(6)); sinLb2 = sin(Eb2(6));
    cosLb1 = cos(Eb1(6)); cosLb2 = cos(Eb2(6));

else
    
    % Initialize estimates for expansion-center states and Jacobians
    % using the mean states and Jacobians provided as input
    
    [Eb1,C2E] = convert_cartesian_to_ephoffset(t,xb1,xb1,eph1,1);
    M3   = [C2E.Xhat_tau  ; C2E.Yhat_tau  ; C2E.Zhat_tau];
    Jb1 = [M3 Z3x3; Z3x3 M3]; % dE/dX
    Qb1 = Jb1 * Pb1 * Jb1';
    
    [Eb2,C2E] = convert_cartesian_to_ephoffset(t,xb2,xb2,eph2,1);
    M3   = [C2E.Xhat_tau  ; C2E.Yhat_tau  ; C2E.Zhat_tau];
    Jb2 = [M3 Z3x3; Z3x3 M3]; % dE/dX
    Qb2 = Jb2 * Pb2 * Jb2';
    
    % Convert dE/dX Jacobians into dX/dE Jacobians
    Jb1 = Jb1'; Jb2 = Jb2';
    
end

% Initialize estimates for expansion-center states and Jacobians
% using the mean states and Jacobians provided as input
xs1 = xb1; Js1 = Jb1; Es1 = Eb1;
xs2 = xb2; Js2 = Jb2; Es2 = Eb2;

% Initialize iteration and convergence variables
iterating = true; iteration = 0; converged = false; failure = 0;
mup_old = NaN; mup_old2 = NaN;

% Iterate to find the peak overlap position and related quantities

while iterating

    if params.verbose
        disp(['iter = ' num2str(iteration,'%03i')]);
        if iteration == 0
            disp([' Eb1 = ' num2str(Es1')]);
            disp([' Eb2 = ' num2str(Es2')]);
        end
    end

    % Difference from mean elements, and offset states for
    % current iteration
    if iteration == 0
        dE1 = zeros(size(Eb1));
        dE2 = zeros(size(Eb2));
        xu1 = xs1;
        xu2 = xs2;
    else
        dE1 = Eb1-Es1;
        dE2 = Eb2-Es2;
        if UseEquinoctialElements
            dE1(6) = asin(sinLb1*cos(Es1(6))-cosLb1*sin(Es1(6)));
            dE2(6) = asin(sinLb2*cos(Es2(6))-cosLb2*sin(Es2(6)));
        end
        xu1 = xs1+Js1*dE1;
        xu2 = xs2+Js2*dE2; 
    end

    % Extract pos & vel vectors from offset states
    ru1 = xu1(1:3); vu1 = xu1(4:6);
    ru2 = xu2(1:3); vu2 = xu2(4:6);

    % Calculate pos/vel covariances, enforcing symmetry to account
    % for round-off errors
    % Slower version
    % aux.Ps1 = cov_make_symmetric( Js1 * Qb1 * Js1' );
    % aux.Ps2 = cov_make_symmetric( Js2 * Qb2 * Js2' );
    % Faster version
    aux.Ps1 = Js1 * Qb1 * Js1'; aux.Ps1 = (aux.Ps1+aux.Ps1')/2;
    aux.Ps2 = Js2 * Qb2 * Js2'; aux.Ps2 = (aux.Ps2+aux.Ps2')/2;

    % Decompose pos/vel covariances into 3x3 submatrices
    As1 = aux.Ps1(1:3,1:3); Bs1 = aux.Ps1(4:6,1:3); % Cs1 = aux.Ps1(4:6,4:6);
    As2 = aux.Ps2(1:3,1:3); Bs2 = aux.Ps2(4:6,1:3); % Cs2 = aux.Ps2(4:6,4:6);

    % As1pAs2 = As1+As2;

    % Calculate inverses of position covariances
    % Slower version
    % [~,~,~,~,CL1,~,As1inv] = CovRemEigValClip(As1,Lclip);
    % [~,~,~,~,CL2,~,As2inv] = CovRemEigValClip(As2,Lclip);
    % if params.verbose
    %     if CL1; disp(' As1 clipped'); end
    %     if CL2; disp(' As2 clipped'); end
    % end
    % Faster version
    [Veig,Leig] = eig(As1); Leig = diag(Leig); 
    Leig(Leig < Lclip) = Lclip;
    As1inv = Veig * diag(1./Leig) * Veig';
    [Veig,Leig] = eig(As2); Leig = diag(Leig); 
    Leig(Leig < Lclip) = Lclip;
    As2inv = Veig * diag(1./Leig) * Veig';

    % Calculate the covariance and peak position of the
    % overlap distribution
    Sigpinv = As1inv + As2inv;
    Sigp = Sigpinv \ I3x3;
    mup = Sigp * ( As1inv*ru1 + As2inv*ru2 );

    % No SIGMAp remediation
    % SigpReminv = Sigpinv \ I3x3;

    % Remediate convergence by convolving SIGMAp matrices with
    % sphere up to radius HBR
    SigpRem = Sigp + SigpRem0;
    SigpReminv = SigpRem \ I3x3;

    % Iterative processing

    if params.maxiter <= 1

        % Iterations forcibly discontinued after one iteration.
        % This yields the original Coppola (2012) 3D-Pc result
        iterating = false;
        converged = true;

        % C12a estimates for vu1-prime and vu2-prime
        vu1p = vu1;
        vu2p = vu2;

    else

        % Do mu-point averaging for iterations past max limit to
        % accelerate slow convergence cases
        if iteration > avgiter
            mup = 0.5*(mup+mup_old);
        end

        % New estimates for vu1-prime and vu2-prime
        vu1p = vu1 + Bs1*As1inv*(mup-ru1);
        vu2p = vu2 + Bs2*As2inv*(mup-ru2);

        % Calculate the energy of the (mup,vu1p) and (mup,vu2p) states
        mupmag = norm(mup);
        Energy0 = -GM/mupmag;
        Energy1 = vu1p'*vu1p/2 + Energy0;
        Energy2 = vu2p'*vu2p/2 + Energy0;

        if params.verbose
            disp([' vu1p = ' num2str(vu1p')]);
            disp([' vu2p = ' num2str(vu2p')]);
            rs1 = xs1(1:3); d1 = (mup-rs1); Msq1 = d1'*As1inv*d1;
            rs2 = xs2(1:3); d2 = (mup-rs2); Msq2 = d2'*As2inv*d2;
            disp([' Msq1,2 = ' num2str(Msq1) ' ' num2str(Msq2)]);
            disp([' Energy1,2 = ' num2str(Energy1) ' ' num2str(Energy2)]);
        end


        % Mark as unconverged if an unbound orbit was
        % encountered during the iteration process

        if UseEquinoctialElements && max(Energy1,Energy2) >= 0

            % Convergence failure due to unbound primary or secondary orbit
            iterating = false;
            converged = false;
            failure = 10*(Energy1 >= 0)+(Energy2 >= 0);

        else

            % New estimates for expansion-center cartesian states
            xs1 = [mup; vu1p];
            xs2 = [mup; vu2p];

            % Iteration and convergence processing

            if iteration > 0

                % Check for convergence using Maha.distance test, imposing
                % maximum iterations

                dmup = mup_old-mup;
                dMD2 = dmup' * SigpReminv * dmup;

                if dMD2 <= params.MD2tol(1)
                    iterating = false;
                    converged = true;
                    dMD2osc   = NaN;
                elseif iteration >= params.maxiter
                    iterating = false;
                    converged = false;
                    dMD2osc   = NaN;
                elseif iteration > osciter
                    % Check for back and forth oscillating convergence
                    dmuposc = mup_old2-mup;
                    dMD2osc = dmuposc' * SigpReminv * dmuposc;
                    if dMD2osc <= params.MD2tol(1)
                        iterating = false;
                        converged = 0.5;
                    else
                        if iteration > avgiter
                            MD2cut = params.MD2tol(2);
                            omfrc = 0;
                        elseif iteration <= acciter
                            MD2cut = params.MD2tol(1);
                            omfrc = 1;
                        else
                            frc = (iteration-acciter)/(avgiter-acciter);
                            omfrc = 1-frc;
                            MD2cut = exp( ...
                                omfrc*log(params.MD2tol(1)) + ...
                                frc  *log(params.MD2tol(2)) );
                        end
                        if dMD2 <= MD2cut
                            iterating = false;
                            converged = 0.1+omfrc/4;
                        end
                    end
                end

                if params.verbose
                    disp([ ...
                        ' rsdiff = ' num2str(norm(mup_old-mup)) ...
                        ' vs1dif = ' num2str(norm(vu1p_old-vu1p)) ...
                        ' vs2dif = ' num2str(norm(vu2p_old-vu2p)) ...
                        ' dMD2 = ' num2str(dMD2) ...
                        ' dMD2osc = ' num2str(dMD2osc)]);
                end

            end

            if UseEquinoctialElements
            
                % Calculate equinoctial elements for primary and
                % secondary at the current iteration's estimate for the
                % expansion-center cartesian states (xs1,xs2).
                [a1s,n1s,af1s,ag1s,chi1s,psi1s,lM1s] = ...
                    convert_cartesian_to_equinoctial(xs1(1:3),xs1(4:6),[],[],verbose);
                [a2s,n2s,af2s,ag2s,chi2s,psi2s,lM2s] = ...
                    convert_cartesian_to_equinoctial(xs2(1:3),xs2(4:6),[],[],verbose);
                
                % Check if any equnoctial orbital elements of the
                % primary and secondary POP states are bad,
                % indicating an unconverged orbit
                if isempty(n1s) || isnan(a1s)
                    bad1s = true;
                else
                    bad1s = false;
                end
                if isempty(n2s) || isnan(a2s)
                    bad2s = true;
                else
                    bad2s = false;
                end
                
                % Check for failures due to unconverged or
                % unbound equinoctial orbit(s)
                if bad1s || bad2s
                    
                    % Convergence failure due to unconverged equinoctial
                    % orbit(s)
                    iterating = false;
                    converged = false;
                    failure = 1e3*bad1s + 1e2*bad2s; % 1100, 1000, or 0100
                    
                else
                    
                    % Check for unbound orbits
                    esq1s = af1s^2+ag1s^2;
                    unbound1s = (a1s <= 0) | esq1s >= 1;
                    esq2s = af2s^2+ag2s^2;
                    unbound2s = (a2s <= 0) | esq2s >= 1;

                    if unbound1s || unbound2s

                        % Convergence failure due to a <= 0 unbound orbit
                        iterating = false;
                        converged = false;
                        failure = 10*unbound1s + unbound2s; % 11, 10, or 01

                    else

                        % Epoch mean longitudes at the initial times
                        lM1s = mod(lM1s,twopi);
                        lM2s = mod(lM2s,twopi);

                        % Equinoctial states, as required for the
                        % next iteration, or for the aux. output
                        Es1 = [n1s;af1s;ag1s;chi1s;psi1s;lM1s];
                        Es2 = [n2s;af2s;ag2s;chi2s;psi2s;lM2s];

                        if params.verbose
                            disp([' Es1 = ' num2str(Es1')]);
                            disp([' Es2 = ' num2str(Es2')]);
                        end

                        % Epoch equinoctial Jacobians as required for the
                        % next iteration, or for the aux. output
                        Js1 = jacobian_E0_to_Xt(0,Es1);
                        Js2 = jacobian_E0_to_Xt(0,Es2);

                        % Save results from this iteration to use in the next
                        if iterating
                            % Comparison variables
                            mup_old2 = mup_old;
                            mup_old  = mup;
                            vu1p_old = vu1p;
                            vu2p_old = vu2p;
                            % MD2_old  = MD2;
                            % Increment iteration counter
                            iteration = iteration+1;
                        end
                        
                    end

                end
                
            else
                
                % Eph. offset states and associated Jacobians, as
                % required for the next iteration or for aux output
                [Es1,AuxEs1] = convert_cartesian_to_ephoffset(t,xs1,xb1,eph1,1);
                [Es2,AuxEs2] = convert_cartesian_to_ephoffset(t,xs2,xb2,eph2,1);
                
                if UseFasterJacobianAlg
                    Js1 = JacobianEphoffset2CartesianFast(t,Es1,AuxEs1,xb1,eph1);
                    Js2 = JacobianEphoffset2CartesianFast(t,Es2,AuxEs2,xb2,eph2);
                else
                    Js1 = JacobianEphoffset2CartesianFull(t,Es1,xb1,eph1);
                    Js2 = JacobianEphoffset2CartesianFull(t,Es2,xb2,eph2);
                end
                
                if params.verbose
                    disp([' Es1 = ' num2str(Es1')]);
                    disp([' Es2 = ' num2str(Es2')]);
                end
                
                % Save results from this iteration to use in the next
                if iterating
                    % Comparison variables
                    mup_old2 = mup_old;
                    mup_old  = mup;
                    vu1p_old = vu1p;
                    vu2p_old = vu2p;
                    % MD2_old  = MD2;
                    % Increment iteration counter
                    iteration = iteration+1;
                end

            end

        end

    end

end

% Assemble the output quantities

conv = converged;
rpk  = mup;
v1pk = vu1p;
v2pk = vu2p;

if Nargout > 4

    aux.converged = converged;
    aux.iteration = iteration;
    aux.failure   = failure;

    if converged

        % Expansion-center pos/vel states and Jacobians
        aux.xs1 = xs1; aux.Js1 = Js1; aux.Es1 = Es1;
        aux.xs2 = xs2; aux.Js2 = Js2; aux.Es2 = Es2;

        % Peak-overlap volume covariance
        aux.Sigp = Sigp; aux.Sigpinv = Sigpinv;
        aux.SigpRem = SigpRem; aux.SigpReminv = SigpReminv;

        if params.maxiter <= 1

            % Original Coppola (2012) result
            aux.xu1 = xb1; aux.dE1 = zeros(6,1);
            aux.xu2 = xb2; aux.dE2 = zeros(6,1);

        else

            % Difference from mean epoch eq. elements
            dE1 = Eb1-Es1;
            dE2 = Eb2-Es2;
            if UseEquinoctialElements
                dE1(6) = asin(sinLb1*cos(Es1(6))-cosLb1*sin(Es1(6)));            
                dE2(6) = asin(sinLb2*cos(Es2(6))-cosLb2*sin(Es2(6)));
            end

            % Calculate offset states for final iteration
            aux.xu1 = xs1+Js1*dE1; aux.dE1 = dE1;
            aux.xu2 = xs2+Js2*dE2; aux.dE2 = dE2;

        end

    end

end

% Report results of iterative process
if params.verbose
    disp(['Converged = ' num2str(converged) ' at t = ' num2str(t)]);
    keyboard;
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
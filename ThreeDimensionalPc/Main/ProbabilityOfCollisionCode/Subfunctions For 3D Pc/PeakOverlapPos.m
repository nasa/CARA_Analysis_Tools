function [conv,rpk,v1pk,v2pk,aux] = PeakOverlapPos(t,xb1,Jb1,t01,Eb01,Qb01,xb2,Jb2,t02,Eb02,Qb02,HBR,params)
% PeakOverlapPos - Find the inertial-frame position, rpk, of the point of
% peak overlap of the primary and secondary distributions for time = t,
% along with the mean primary and secondary velocities at that point.
%
% Syntax: [conv,rpk,v1pk,v2pk,aux] = PeakOverlapPos(t,xb1,Jb1,t01,Eb01,Qb01,xb2,Jb2,t02,Eb02,Qb02,HBR,params)
%
% =========================================================================
%
% Description:
%
% Find the inertial-frame position, rpk, of the point of
% peak overlap the primary and secondary distributions for time = t. The
% primary and secondary distributions arise from initial equinoctial 
% MVN distributions defined at times (t01,t02), with mean eq. states and 
% covariances (Eb01,Qb01) and (Eb02,Qb02).
%
% All time units in s and length units in m, unless otherwise noted.
%
% =========================================================================
%
% Inputs:
%
%   t    = Current time (s)
%   xb1  = Mean inertial-frame pos/vel state for primary  [6x1]
%   Jb1  = Jacobian matrix dx/dE evaluated at mean state for primary [6x6]
%   t01  = Initial time for equinoctial PDF for primary [1x1]
%   Eb01 = Initial mean equinoctial element state for primary [6x1]
%          (Note xb1 and Eb01 map to one another through two-body equations
%          of motion).
%   Qb01 = Iniital mean equinoctial covariance for primary [6x6]
%   xb2  = Mean inertial-frame pos/vel state for secondary  [6x1]
%   Jb2  = Jacobian matrix dx/dE evaluated at mean state for secondary [6x6]
%   t02  = Initial time for equinoctial PDF for secondary [1x1]
%   Eb02 = Initial mean equinoctial element state for secondary [6x1]
%          (Note xb2 and Eb02 map to one another through two-body equations
%          of motion).
%   Qb02 = Iniital mean equinoctial covariance for secondary [6x6]
%   HBR  = Combined primary+secondary hard-body radius [1x1]
%   params = Structure of execution parameters (optional - see the
%            comments below for detailed description).
%
% =========================================================================
%
% Outputs:
%
%   conv = Convergence flag (true or false)
%   rpk  = Position of peak overlap of primary & secondary PDFs
%   v1pk = Mean primary velocity at rpk
%   v2pk = Mean secondary velocity at rpk
%   aux  = Auxilliary outputs
%
% =========================================================================
%
% References:
%
%    Doyle T. Hall (2020) "Satellite Collision Rates and Probabilities"
%    in preparation
%
%    Vincent T. Coppola (2012a) "Including Velocity Uncertainty in the
%    Probability of Collision Between Space Objects" AAS 12-247.
%
%    Vincent T. Coppola (2012b) "Evaluating the Short Encounter Assumption
%    of the Probability of Collision Formula" AAS 12-248.
%
%    Hereafter these will be referred to as "H20", "C12a" and "C12b".
%
% =========================================================================
%
% Examples/Validation Cases: None
%
% Other m-files required:
%  cov_make_symmetric.m
%  CovRemEigValClip.m
%  convert_cartesian_to_equinoctial.m
%  jacobian_E0_to_Xt.m
%
% Subfunctions: None
%
% MAT-files required: None
%
% Initial version: Jan 2020
%
% ----------------- BEGIN CODE -----------------

    % Initializations and defaults
    Nargin = nargin; Nargout = nargout;
    
    % Set up default parameters
    if Nargin < 13; params = []; end

    % Tolerance for convergence using relative Maha. distance squared 
    if ~isfield(params,'MD2tol') || isempty(params.MD2tol)
        params.MD2tol = 1e-6;
    end
    
    % Maximum number of iterations to perform
    if ~isfield(params,'maxiter') || isempty(params.maxiter)
        params.maxiter = 100;
    end
    avgiter = min(100,round(params.maxiter/2));

    % Eigenvalue clipping factor
    if ~isfield(params,'Fclip') || isempty(params.Fclip)
        params.Fclip = 1e-4;
    end
    Lclip = (HBR*params.Fclip)^2;
    
    if ~isfield(params,'GM'); params.GM = []; end;
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
    
    % Sin and Cos values for initial mean longitudes
    sinLb01 = sin(Eb01(6)); sinLb02 = sin(Eb02(6));
    cosLb01 = cos(Eb01(6)); cosLb02 = cos(Eb02(6));
    
    % Offsets of current time from the initial times of equinoctial states
    dt01 = t-t01; dt02 = t-t02;
    
    % Initialize estimates for expansion-center states and Jacobians
    % using the mean states and Jacobians provided as input
    xs1 = xb1; Js1 = Jb1; Es01 = Eb01;
    xs2 = xb2; Js2 = Jb2; Es02 = Eb02;

    % Initialize iteration and convergence variables
    iterating = true; iteration = 0; converged = false; failure = 0;

    % Iterate to find the peak overlap position and related quantities

    while iterating
        
        if params.verbose
            disp(['iter = ' num2str(iteration,'%03i')]);
            if iteration == 0
                disp([' Eb01 = ' num2str(Es01')]);
                disp([' Eb02 = ' num2str(Es02')]);
            end
        end

        % Difference from mean epoch eq. elements, and offset states for
        % current iteration
        if iteration == 0
            dE01 = zeros(size(Eb01));
            dE02 = zeros(size(Eb02));
            xu1 = xs1;
            xu2 = xs2;
        else
            dE01(1:5) = Eb01(1:5)-Es01(1:5);
            dE01(6) = asin(sinLb01*cos(Es01(6))-cosLb01*sin(Es01(6)));
            dE02(1:5) = Eb02(1:5)-Es02(1:5);
            dE02(6) = asin(sinLb02*cos(Es02(6))-cosLb02*sin(Es02(6)));
            xu1 = xs1+Js1*dE01;
            xu2 = xs2+Js2*dE02; 
        end

        % Extract pos & vel vectors from offset states
        ru1 = xu1(1:3); vu1 = xu1(4:6);
        ru2 = xu2(1:3); vu2 = xu2(4:6);

        % Calculate pos/vel covariances, enforcing symmetry to account
        % for round-off errors
        Ps1 = cov_make_symmetric( Js1 * Qb01 * Js1' );
        Ps2 = cov_make_symmetric( Js2 * Qb02 * Js2' );

        % Decompose pos/vel covariances into 3x3 submatrices
        As1 = Ps1(1:3,1:3); Bs1 = Ps1(4:6,1:3); % Cs1 = Ps1(4:6,4:6);
        As2 = Ps2(1:3,1:3); Bs2 = Ps2(4:6,1:3); % Cs2 = Ps2(4:6,4:6);

        % Calculate inverses of position covariances
        [~,~,~,~,CL1,~,As1inv] = CovRemEigValClip(As1,Lclip);
        [~,~,~,~,CL2,~,As2inv] = CovRemEigValClip(As2,Lclip);
        if params.verbose
            if CL1; disp(' As1 clipped'); end
            if CL2; disp(' As2 clipped'); end
        end
        
        % Calculate the covariance and peak position of the
        % overlap distribution
        Sigpinv = As1inv + As2inv;
        Sigp = Sigpinv \ I3x3;
        mup = Sigp * ( As1inv*ru1 + As2inv*ru2 );
        
        if params.verbose
            disp([' mup = ' num2str(mup')]);
        end

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

            if (Energy1 >= 0) || (Energy2 >= 0)

                % Convergence failure due to unbound primary or secondary orbit
                iterating = false;
                converged = false;
                failure = 10*(Energy2 >= 0)+(Energy2 >= 0);

            else
                
                % New estimates for expansion-center cartesian states
                xs1 = [mup; vu1p];
                xs2 = [mup; vu2p];

                % % Calculate relative position covariance and inverse
                % As = As1+As2;
                % [~,~,~,~,~,~,Asinv] = CovRemEigValClip(As,Lclip);
                % 
                % % Calculate rel.pos. Maha.distance squared
                % ru = ru2-ru1;
                % MD2 = ru' * Asinv * ru;

                % Iteration and convergence processing

                if iteration > 0

                    % Check for convergence using Maha.distance test, imposing
                    % maximum iterations

                    dmup = mup_old-mup;
                    dMD2 = dmup' * Sigpinv * dmup;
                    % dMD2 = dmup' * Asinv * dmup;

                    if params.verbose
                        disp([ ...
                            ' rsdiff = ' num2str(norm(mup_old-mup)) ...
                            ' vs1dif = ' num2str(norm(vu1p_old-vu1p)) ...
                            ' vs2dif = ' num2str(norm(vu2p_old-vu2p)) ...
                            ' dMD2   = ' num2str(dMD2)]);
                    end

                    if (dMD2 <= params.MD2tol)
                        iterating = false;
                        converged = true;
                    elseif iteration >= params.maxiter
                        iterating = false;
                        converged = false;
                    end

                    % MD2_dif = abs(MD2-MD2_old);
                    % if (MD2_dif <= params.MD2tol)
                    %     iterating = false;
                    %     converged = true;
                    % elseif iteration >= params.maxiter
                    %     iterating = false;
                    %     converged = false;
                    % end
                    % 
                    % if params.verbose
                    %     disp([ ...
                    %         ' rsdiff = ' num2str(norm(mup_old-mup)) ...
                    %         ' vs1dif = ' num2str(norm(vu1p_old-vu1p)) ...
                    %         ' vs2dif = ' num2str(norm(vu2p_old-vu2p)) ...
                    %         ' MD2del = ' num2str(dMD2) ...
                    %         ' MD2dif = ' num2str(MD2_dif)]);
                    % end

                end

                % Calculate equinoctial elements for primary and
                % secondary at the current iteration's estimate for the
                % expansion-center cartesian states (xs1,xs2).
                [a1s,n1s,af1s,ag1s,chi1s,psi1s,lM1s] = ...
                    convert_cartesian_to_equinoctial(xs1(1:3),xs1(4:6),[],[],verbose);
                [a2s,n2s,af2s,ag2s,chi2s,psi2s,lM2s] = ...
                    convert_cartesian_to_equinoctial(xs2(1:3),xs2(4:6),[],[],verbose);

                % Mark as unconverged if an unbound orbit was
                % encountered during the iteration process

                if (a1s <= 0) || (a2s <= 0)

                    % Convergence failure due to unbound primary or secondary orbit
                    iterating = false;
                    converged = false;
                    failure = 10*(a1s <= 0)+(a2s <= 0);

                else

                    % Epoch mean longitudes at the initial times
                    lM10s = mod(lM1s-n1s*dt01,twopi);
                    lM20s = mod(lM2s-n2s*dt02,twopi);

                    % Epoch equinoctial states and associated Jacobians, as
                    % required for the next iteration or for aux output
                    Es01 = [n1s;af1s;ag1s;chi1s;psi1s;lM10s];
                    Es02 = [n2s;af2s;ag2s;chi2s;psi2s;lM20s];            
                    % [~,~,~,Js1] = jacobian_E0_to_Xt(dt01,Es01,[],[]);
                    % [~,~,~,Js2] = jacobian_E0_to_Xt(dt02,Es02,[],[]);
                    Js1 = jacobian_E0_to_Xt(dt01,Es01);
                    Js2 = jacobian_E0_to_Xt(dt02,Es02);

                    if params.verbose
                        disp([' Es01 = ' num2str(Es01')]);
                        disp([' Es02 = ' num2str(Es02')]);
                    end

                    % Save results from this iteration to use in the next

                    if iterating
                        % Comparison variables
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
            aux.xs1 = xs1; aux.Js1 = Js1; aux.Es01 = Es01;
            aux.xs2 = xs2; aux.Js2 = Js2; aux.Es02 = Es02;
            
            % Peak-overlap volume covariance
            aux.Sigp = Sigp; aux.Sigpinv = Sigpinv;
            
            if params.maxiter <= 1
                
                % Original Coppola (2012) result
                aux.xu1 = xb1; aux.dE01 = zeros(6,1);
                aux.xu2 = xb2; aux.dE02 = zeros(6,1);
                
            else
                
                
                % Difference from mean epoch eq. elements
                dE01(1:5) = Eb01(1:5)-Es01(1:5);
                dE01(6) = asin(sinLb01*cos(Es01(6))-cosLb01*sin(Es01(6)));
                dE02(1:5) = Eb02(1:5)-Es02(1:5);
                dE02(6) = asin(sinLb02*cos(Es02(6))-cosLb02*sin(Es02(6)));

                % Calculate offset states for final iteration
                aux.xu1 = xs1+Js1*dE01; aux.dE01 = dE01;
                aux.xu2 = xs2+Js2*dE02; aux.dE02 = dE02;
                
            end

        end
        
    end

    % Report results of iterative process
    if params.verbose
        disp(['Converged = ' num2str(converged) ' at t = ' num2str(t)]);
    end

    % if ~converged
    %     keyboard;
    % end

    return
    
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D.Hall         | 2020-JAN-10 | Initial Development
% D.Hall         | 2020-JUN-16 | Modified to allow maxiter = 1 to calculate
%                                the original Coppola (2012) formulation.
function [PX,PXA,PXB,PE,PEA,PEB] = TBSTMCovInterp(T,X,T1,X1,PX1,T2,X2,PX2,method,mu)
%
% TBSTMEstimCov - Use 2-body state transition matrices (TBSTMs) to estimate
% an interpolated cartesian/inertial-frame (pos,vel) covariance, PX,
% at a time = T, from the bracketing (pos,vel) states and covariances
% (T1,X1,PX1) and (T2,X2,PX2).
%
% Syntax: [PX,PXA,PXB,PE,PEA,PEB] = TBSTMCovInterp(T,X,T1,X1,PX1,T2,X2,PX2,method,mu)
%
% [PX,PXA,PXB,PE,PEA,PEB] = TBSTMEstimCov(T,X,T1,X1,PX1,T2,X2,PX2,method,mu)
%
% Inputs:
%
%   T           = Time for output interpolated covariance (s) [1x1]
%   X           = Cartesian (pos,vel) state at time = T (km, km/s) [6x1]
%                 NOTE: Inputing a pre-interpolated state here is optimal.
%                       However, if X is input empty, then this function
%                       will interpolate a rough estimate for X for itself.
%   (T1,X1,PX1) = First bracketing time, state, & covariance
%                  T1:  (s) [1x1] 
%                  X1:  (km, km/s) [6x1]
%                  PX1: (km^2, km^2/s, (km/s)^2) [6x6]
%   (T2,X2,PX2) = Second bracketing time, state, & covariance
%                  T2:  (s) [1x1] 
%                  X2:  (km, km/s) [6x1]
%                  PX2: (km^2, km^2/s, (km/s)^2) [6x6]
%   method      = Two-body state transistion matrix method:
%                  'cartesian'   = Cartesian STM method (Shepperd, 1985)
%                  'equinoctial' = Equinoctial STM method
%                  DEFAULT = 'cartesian'
%   mu          = Gravitational constant
%                  DEFAULT = EGM-96 value [km^3/s^2]
%
% Outputs:
%
%   PX          = Interpolated cartesian covariance
%                  (km^2, km^2/s, (km/s)^2) [6x6]
%   PXA         = Forward-propagated cartesian covariance  (T1 -> T)
%                  (km^2, km^2/s, (km/s)^2) [6x6]
%   PXB         = Backward-propagated cartesian covariance (T2 -> T)
%                  (km^2, km^2/s, (km/s)^2) [6x6]
%
% Optional outputs (request only as needed to save CPU time):
%
%   PE          = Interpolated equinoctial covariance [6x6]
%                  Using Vallado and Alfano (2015) ordering
%                  [N,Af,Ag,Chi,Psi,L]
%   PEA         = Forward-propagated equinoctial covariance [6x6]
%   PEB         = Backward-propagated equinoctial covariance [6x6]
%
% Other m-files required:
%     k2b_state_transition.m
%     convert_cartesian_to_equinoctial.m
%     convert_equinoctial_to_cartesian.m
%     jacobian_equinoctial_to_cartesian
%   
% Subfunctions: None
% MAT-files required: None
%
% See also:
%    D. T. Hall, "Orbital State Covariance Interpolation Using State 
%    Transition Matrices" Omitron Tech. Document, 2019 Dec 16.
%
% Initial version: Aug 2019 ; Last Revision: Feb 2023
%
% ----------------- BEGIN CODE -----------------

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p, '../OrbitTransformations')); addpath(s.path);
    pathsAdded = true;
end

% Initializations

Nargin = nargin; Nargout = nargout;

if Nargin < 10 || isempty(mu)
    % Earth gravitational constant (EGM-96) [km^3/s^2]
    mu  = 3.986004418e5;
end

if (Nargin < 9) || isempty(method)
    % Default to equinoctial element 2-body motion STM formulation
    method = 'equinoctial';
    % Cartesian STM formulation of Shepperd (1984) is slightly faster but
    % occaisonally less accurate
    % method = 'cartesian';
end

% Linear interpolation weights

tau = (T-T1)/(T2-T1); omtau = 1-tau;

% Use specified STM estimation method

method = lower(method);

switch method

    case 'cartesian'
        
        % Calculate the two-body STMs using the Shepperd (1985) method
        
        [~,~,XA,STMA] = k2b_state_transition(X1,T-T1,mu,[]);
        [~,~,XB,STMB] = k2b_state_transition(X2,T-T2,mu,[]);
        
        % Calculate the forward- and backward-propagated covariances

        PXA = STMA * PX1 * STMA';
        PXB = STMB * PX2 * STMB';
        
        % Weighted covariance interpolant
        
        PX = omtau*PXA+tau*PXB;
        
        % Only calculate equinoctial covariances if required
        
        if Nargout > 3
            
            % Weighted state interpolant (if required)

            if isempty(X)
                X = omtau*XA+tau*XB;
            end

            % Equinoctial state at interpolated time
            
            [~,n,af,ag,chi,psi,lam,~] = ...
                convert_cartesian_to_equinoctial(X(1:3),X(4:6),1,mu);
            E = [n; af; ag; chi; psi; lam];
            
            % Jacobian (dX/dE) at interpolated time

            Jetoc = jacobian_equinoctial_to_cartesian(E,X,1,mu);
            Jctoe = Jetoc\eye(6,6);
            
            % Equinoctial covariances at interpolated time
            
            PE = Jctoe * PX * Jctoe';
            
            if Nargout > 4
                PEA = Jctoe * PXA * Jctoe';
                PEB = Jctoe * PXB * Jctoe';
            end

        end
        
    case 'equinoctial'
        
        % Convert to equinoctial state representation

        % Equinoctial states

        [~,n1,af1,ag1,chi1,psi1,lam1,~] = ...
            convert_cartesian_to_equinoctial(X1(1:3),X1(4:6),1,mu);
        E1 = [n1; af1; ag1; chi1; psi1; lam1];

        [~,n2,af2,ag2,chi2,psi2,lam2,~] = ...
            convert_cartesian_to_equinoctial(X2(1:3),X2(4:6),1,mu);
        E2 = [n2; af2; ag2; chi2; psi2; lam2];

        % Jacobians

        I6x6 = eye(6,6);
        
        Jetoc1 = jacobian_equinoctial_to_cartesian(E1,X1,1,mu);
        Jctoe1 = Jetoc1\I6x6;
        % Jctoe1 = inv(Jetoc1);

        Jetoc2 = jacobian_equinoctial_to_cartesian(E2,X2,1,mu);
        Jctoe2 = Jetoc2\I6x6;
        % Jctoe2 = inv(Jetoc2);

        % Covariances

        PE1 = Jctoe1 * PX1 * Jctoe1';  
        PE2 = Jctoe2 * PX2 * Jctoe2';
        
        % Equinoctial 2-body motion STMs only have one off-diagonal element
        
        STM1 = I6x6; STM1(6,1) = T-T1;
        STM2 = I6x6; STM2(6,1) = T-T2;
        
        % Calculate the forward- and backward-propagated covariances

        PEA = STM1 * PE1 * STM1';
        PEB = STM2 * PE2 * STM2';
        
        % Weighted equinoctial covariance interpolant
        
        PE = omtau*PEA+tau*PEB;
        
        if isempty(X)
            % Weighted state interpolant (if required)
            E = omtau*E1+tau*E2;
            % Interpolate lambda in sine and cosine space
            lamA = lam1+n1*(T-T1); clamA = cos(lamA); slamA = sin(lamA);
            lamB = lam2+n2*(T-T2); clamB = cos(lamB); slamB = sin(lamB);
            E(6) = atan2(omtau*slamA+tau*slamB,omtau*clamA+tau*clamB);
            [rvec,vvec] = convert_equinoctial_to_cartesian( ...
                E(1),E(2),E(3),E(4),E(5),E(6),0,1,mu);
            X = [rvec; vvec];
        else
            % Equinoctial state at interpolated time
            [~,n,af,ag,chi,psi,lam,~] = ...
                convert_cartesian_to_equinoctial(X(1:3),X(4:6),1,mu);
            E = [n; af; ag; chi; psi; lam];
            if isempty(E)
                error('Invalid cartesian to equinoctial conversion');
            end
        end

        % Jacobian (dX/dE) at interpolated time

        Jetoc = jacobian_equinoctial_to_cartesian(E,X,1,mu);
        
        % Cartesian covariances at interpolated time
        
        PX = Jetoc * PE * Jetoc';
        
        if Nargout > 1
            PXA = Jetoc * PEA * Jetoc';
            PXB = Jetoc * PEB * Jetoc';
        end

    otherwise
        
        error('Invalid STM estimation method');
        
end

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
% D.Hall         | 2019-AUG-27 | Initial Development
% D.Hall         | 2019-DEC-16 | Completed testing and documentation
% L. Baars       | 2022-OCT-03 | Fixed pathing for SDK restructuring
% L. Baars       | 2023-FEB-27 | Fixed relative pathing issue in addpath
%                                calls.
%


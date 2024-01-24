function [X,P,E,J,K,Q,QRemStat,QRaw,QRem,CRem] = EquinoctialMatrices(r,v,C,RemEqCov)

% Calculate equinoctial state, covariance and Jacobian matrices

% Cartestian state vector and covariance matrix
X = [r; v]/1e3; % Pos/vel state vec. in km units
P = C/1e6;      % Pos/vel covariance in km units

% Calculate the equinoctial state vector and associated matrices
try
    % Calculate equinoctial elements from the pos-vel state
    [~,n,af,ag,chi,psi,lM,~] = ...
        convert_cartesian_to_equinoctial(X(1:3),X(4:6));
    lM = mod(lM,2*pi);
    % Equinoctial state vector, E
    E = [n,af,ag,chi,psi,lM]';
    % Jacobian matrix, J = dX/dE
    J = jacobian_equinoctial_to_cartesian(E,X);
    % Inverse of Jacobian matrix, K = dE/dX
    K = J \ eye(6,6);
    % Equinoctial state covariance, Q = K * P * K'
    Q = cov_make_symmetric(K * P * K');
    % Save raw covariance
    QRaw = Q;
    % Remediate eq. covariance, if required
    if RemEqCov
        % Calc Q remediation status and remediated Q matrix
        [~,~,~,~,QRemStat,~,~,QRem] = ...
            CovRemEigValClip(Q);
        if QRemStat
            Q = cov_make_symmetric(QRem);
            P = cov_make_symmetric(J * Q * J');
            CRem = 1e6*P;
        else
            CRem = C;
        end        
    else
        % Calc Q remediation status only
        [~,~,~,~,QRemStat] = CovRemEigValClip(Q);
        QRem = Q;
        CRem = C;
    end
catch
    warning('EquinoctialMatrices:CalculationFailure','Equinoctial state/covariance/Jacobian calculation failure');
    E        = NaN(6,1);
    J        = NaN(6,6);
    K        = NaN(6,6);
    Q        = NaN(6,6);
    QRemStat = NaN(1,1);    
    QRaw     = NaN(6,6);
    QRem     = NaN(6,6);
    if RemEqCov
        CRem = NaN(6,6);
    else
        CRem = C;
    end
end

return
end

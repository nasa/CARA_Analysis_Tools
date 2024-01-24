function [EA] = Mean2EccAnomaly(M,e)

    sz = numel(M);
    EA = zeros(sz,1);

    id = -pi<M & M<0 | M>pi;
    EA(id)  = M(id) - e(id);
    % EA(~id) = M(~id) + e(~id);
    id = ~id;
    EA(id)  = M(id) + e(id);

    % if (-pi<M && M<0 || M>pi)
    %     EA = M-e;
    % else
    %     EA = M+e;
    % end
    
    tol   = 1e-13 * ones(sz,1);
    delta = ones(sz,1);
    
    EAold = EA;
    EAnew = EA;
    
    while(delta > tol)
        
        EAnew = EAold + ((M-EAold+e.*sin(EAold))./(1-e.*cos(EAold)));
        
        delta = abs(EAnew - EAold);
        
        EAold = EAnew;
    
    end
    
    EA = EAnew;
    
end
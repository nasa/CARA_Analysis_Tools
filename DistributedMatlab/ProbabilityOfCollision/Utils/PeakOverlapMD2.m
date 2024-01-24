function [MD2,Xu,Ps,Asdet,Asinv,POPconverged] = PeakOverlapMD2(t,t10,Eb10,Qb10,t20,Eb20,Qb20,HBR,EMD2,params)

% Calculate the effective Maha. distance using the peak overlap point.
%
% Actually, calculate MD^2 + log(det(As)), to account for slight variations
% in the determinant of the As matrix, which shows up in the denominator.

% Find the peak overlap position (POP) at the nominal TCA

% Calc pos/vel mean states and associated Jacobians at specified time
[Jb1,xb1] = jacobian_E0_to_Xt(t,Eb10);
[Jb2,xb2] = jacobian_E0_to_Xt(t,Eb20);

% Calc peak overlap position
% params.verbose = true;
[POPconverged,~,~,~,POP] = PeakOverlapPos(t, ...
    xb1,Jb1,t10,Eb10,Qb10, ...
    xb2,Jb2,t20,Eb20,Qb20, ...
    HBR,params);

if POPconverged > 0
    
    % Relative POP-corrected distance
    Xu = POP.xu2-POP.xu1; ru = Xu(1:3);
    
    % Construct the POP-corrected joint pos/vel covariance
    if params.XCprocessing
        % DCP 6x1 sensitivity vectors for cartesian pos/vel
        % state at the current eph time - see Hall (2021) eq 57
        GCp = POP.Js1*params.GEp;
        GCs = POP.Js2*params.GEs;
        Ps  = POP.Ps1 + POP.Ps2 - params.sigpXsigs * (GCs*GCp'+GCp*GCs');
    else
        Ps  = POP.Ps1 + POP.Ps2;
    end
    
    % Calculate the inverse of A, remediating NPDs with eigenvalue clipping
    As = Ps(1:3,1:3);
    % Lclip = (params.Fclip*HBR)^2;
    Lclip = 0; % Required for one Obj 48901 DCP case
    % Slower version
    % [~,~,~,~,~,Asdet,Asinv] = CovRemEigValClip(As,Lclip);
    % Faster version
    [Veig,Leig] = eig(As); Leig = diag(Leig); Leig(Leig < Lclip) = Lclip;
    Asdet = Leig(1) * Leig(2) * Leig(3);
    Asinv = Veig * diag(1./Leig) * Veig';
    
    % Calculate effective Maha distance
    MD2 = ru' * Asinv * ru;
    
    % Calculate the "effective" or modified Maha distance,
    % by adding ln(|As|), or ln(|As|) - 2*ln(|vu])
    if EMD2 == 1
        MD2 = MD2 + log(Asdet);
    elseif EMD2 == 2
        MD2 = MD2 + log(Asdet) - 2*log(norm(Xu(4:6)));
    end
    
else
    
    % Return null values for no POP convergence
    MD2 = NaN; Xu = []; Ps = []; Asdet = []; Asinv = [];
    
end

return
end
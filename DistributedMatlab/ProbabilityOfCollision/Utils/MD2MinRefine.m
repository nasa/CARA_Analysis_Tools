function [out] = MD2MinRefine(tmd,tsg,dtsg,ttol,itermax,findmin,Eb10,Qb10,Eb20,Qb20,HBR,POPPAR)
% MD2MinRefine - Iteratively refine the paramters of the min of effective
%                MD2 curve
%
% Syntax: [out] = MD2MinRefine(tmd,tsg,dtsg,ttol,itermax,findmin,Eb10,Qb10,Eb20,Qb20,HBR,POPPAR)
%
% =========================================================================
%
% Copyright (c) 2022-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% DESCRIPTION:
%
%   Refines key parameters to the minimum of the effective Mahalanobis
%   distance squared (MD2) curve. The key parameters are:
%     - The time of the minimum of the MD2 curve
%     - The value of the MD2 curve at the minimum
%     - The numerically derived first derivative of MD2 at the minimum time
%     - The numerically derived second derivative of MD2 at the minimum
%       time
%     - The 1-sigma time width of the MD2 curve at the minimum time
%     - The inertial-frame relative pos/vel state at the minimum time
%     - The inertial-frame combined covariance at the minimum time
%     - The determinant of the position 3x3 of the combined covariance at
%       the minimum time
%     - The inverse of the position 3x3 of the combined covariance at the
%       minimum time
%
%   The function can either calculate these values at the minimum time
%   passed in, or it can search the n-sigma time span for a minimum of the
%   MD2 curve (starting at the "tmd" time passed in) and calculate the
%   values from that new minimum.
%
% =========================================================================
%
% INPUT:
%
%   tmd = Time, in seconds from TCA, of the minimum of the Mahalanobis
%         distance squared parabola
%   tsg = The 1-sigma time width for numerically computing the output
%         values. It is also defines the search interval when findmin is
%         set to true.
%   dtsg = Sigma multiplier applied to the tsg parameter.
%   ttol = Convergence tolerance (in seconds)
%   itermax = Maximum number of iterations performed by the function.
%   findmin = When set to true, indicates that the function should search
%             for the minimum of the MD2 curve. When false, assume that tmd
%             is already the minimum of the MD2 curve.
%   Eb10 = Equinoctial state of the primary object at TCA. [6x1]
%   Qb10 = Equinoctial state covariance of the primary object at TCA. [6x6]
%   Eb20 = Equinoctial state of the secondary object at TCA. [6x1]
%   Qb20 = Equinoctial state covariance of the secondary object at TCA. [6x6]
%   HBR = Combined HBR of the two objects (in m).
%   POPPAR = Peak Overlap (PeakOveralapMD2.m) input parameters.
%
% =========================================================================
%
% OUTPUT:
%
%   out = Structure containing results and other pertinent parameters
%         calculated by the function. Values include:
%     tmd - Time (in seconds from TCA) of the minimum of the MD2 curve
%     ydot - First derivative of the MD2 curve at tmd
%     ydotdot - Second derivative of the MD2 curve at tmd
%     tsg - 1-sigma time width of the MD2 curve at tmd
%     iter - Number of iterations performed by the function
%     MD2converged - Indicates of MD2MinRefine converged
%     POPconverged - Indicates if the calls to PeakOverlapMD2 converged
%     ymd - Value of PeakOverlapMD2 at tmd
%     delt - Delta time for numerical calculation
%     tlo - Time (in sec relative to TCA) of the lower bracketing point
%     thi - Time (in sec relative to TCA) of the upper bracketing point
%     ylo - Value of PeakOverlapMD2 at tlo
%     Xulo - Inertial-frame relative pos/vel state at tlo [6x1]
%     Pslo - Inertial-frame combined covariance at tlo [6x6]
%     Asdetlo - Determinant of the postion 3x3 of Pslo
%     Asinvlo - Inverse of the position 3x3 of Pslo [3x3]
%     Xu - Inertial-frame relative pos/vel state at tmd [6x1]
%     Ps - Inertial-frame combined covariance at tmd [6x6]
%     Asdet - Determinant of the position 3x3 of Ps
%     Asinv - Inverse of the position 3x3 of Ps [3x3]
%     yhi - Value of PeakOverlapMD2 at thi
%     Xuhi - Inertial-frame relative pos/vel state at thi [6x1]
%     Pshi - Inerital-frame combined covariance at thi [6x6]
%     Asdethi - Determinant of the position 3x3 of Pshi
%     Asinvhi - Inverse of the position 3x3 of Pshi [3x3]
%
% =========================================================================
%
% Initial version: Oct 2022; Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

% Iteratively refine minimum of the MD^2 curve in time
iterating = true; out.iter = 0;
itermin = min(2,itermax); itermed = min(4,itermax);
out.MD2converged = false; out.POPconverged = true;

out.tmd = tmd; out.tsg = tsg; out.ymd = NaN;

while iterating
    
    % Increment iteration counter
    out.iter = out.iter+1;
    
    % Save old value of best MD^2 minimization time
    tmdold = out.tmd; tsgold = out.tsg;
    
    % Calculate 3 points to estimate numerical derivatives
    
    % Point spacing
    out.delt = out.tsg*dtsg; twodelt = 2*out.delt; delt2 = out.delt^2;
    
    % Low and high bracketing points
    out.tlo = out.tmd-out.delt;
    out.thi = out.tmd+out.delt;
    
    % Calculate y = MD^2 for the three points
    [out.ylo,out.Xulo,out.Pslo,out.Asdetlo,out.Asinvlo,~,aulo] = ...
        PeakOverlapMD2(out.tlo, ...
             0,Eb10,Qb10, ...
             0,Eb20,Qb20, ...
             HBR,1,POPPAR);
    if out.iter == 1 || findmin
        [out.ymd,out.Xu,out.Ps,out.Asdet,out.Asinv,~,aumd] = ...
            PeakOverlapMD2(out.tmd, ...
                 0,Eb10,Qb10, ...
                 0,Eb20,Qb20, ...
                 HBR,1,POPPAR);
    end
    [out.yhi,out.Xuhi,out.Pshi,out.Asdethi,out.Asinvhi,~,auhi] = ...
        PeakOverlapMD2(out.thi,...
             0,Eb10,Qb10,    ...
             0,Eb20,Qb20,    ...
             HBR,1,POPPAR);

    if isnan(out.ylo) || isnan(out.ymd) || isnan(out.yhi)
        
        % POP search unconverged for one or more of the three points
        iterating = false; out.POPconverged = false;
        
    else

        % Calculate new value of SigmaT
        out.ydot = (out.yhi-out.ylo)/twodelt;
        out.ydotdot = (out.yhi-2*out.ymd+out.ylo)/delt2;
        out.tsg = sqrt(2/max(0,out.ydotdot));

        % Analyze convergence criteria
        if findmin
            % In this case, estimate the new time that minimizes the
            % effective MD^2 curve before analyzing convergence
            out.tmd = out.tmd - out.ydot/out.ydotdot;
            cnvg = abs(out.tmd-tmdold) < ttol*out.tsg;
        else
            cnvg = abs(out.tsg-tsgold) < ttol*out.tsg;
        end

        % Check for Maha. distance convergence
        if ~cnvg
            if findmin 
                MD2check = out.iter >= itermed;
            else
                MD2check = out.iter >= itermin;
            end
            if MD2check
                MD2conv = 1 - ...
                    [aulo.MD2actual auhi.MD2actual]/aumd.MD2actual;
                if max(abs(MD2conv)) < ttol
                    cnvg = true;
                end
            end
        end
        
        if isinf(out.tsg)
            % Second derivative not positive (i.e., ydotdot <= 0)
            iterating = false; out.MD2converged = false;
        elseif cnvg && out.iter >= itermin
            iterating = false; out.MD2converged = true;
        else
            iterating = out.iter <= itermax;
        end
        
        % % For debugging
        % if findmin
        %     cnvgfact = abs(out.tmd-tmdold)/out.tsg;
        % else
        %     cnvgfact = abs(out.tsg-tsgold)/out.tsg;
        % end
        % disp(['Iter = ' num2str(out.iter,'%02i') ...
        %      ' tmd = ' num2str(out.tmd,'%0.8g') ...
        %      ' ymd = ' num2str(out.ymd,'%0.8g') ...
        %      ' tsg = ' num2str(out.tsg,'%0.8g') ...
        %      ' dmd = ' num2str(abs(out.tmd-tmdold)/out.tsg,'%0.8g') ...
        %      ' dsg = ' num2str(abs(out.tsg-tsgold)/out.tsg,'%0.8g') ... 
        %      ' cnv = ' num2str(cnvgfact)]);
            
    end
    
end

return;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2022-Oct-03 | Initial Development
% D. Hall        | 2023-Mar-01 | Updated method for calculating ydot and
%                                ydotdot. Added convergence check for
%                                ydotdot <= 0.
% D. Hall        | 2025-Jan-02 | Added checks for Mahalanobis distance
%                                convergence.
% L. Baars       | 2025-Aug-06 | Updated documentation for public release.

% =========================================================================
%
% Copyright (c) 2022-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

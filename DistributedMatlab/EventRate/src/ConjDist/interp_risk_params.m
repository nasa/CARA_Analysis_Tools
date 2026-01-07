function out = interp_risk_params(L10Pcmission,L10PremSearch, ...
                                  L10Pcummd,L10Pcumlo,L10Pcumhi, ...
                                  Nmvrmd,Nmvrlo,Nmvrhi)
% interp_risk_params - Interpolate the risk parameters Prem, Nmvr and Pcum
%
% Syntax: out = interp_risk_params(L10Pcmission,L10PremSearch, ...
%                                  L10Pcummd,L10Pcumlo,L10Pcumhi, ...
%                                  Nmvrmd,Nmvrlo,Nmvrhi)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Interpolate the risk parameters Prem, Nmvr and Pcum
%
% =========================================================================
%
% Input:
%
%   L10Pcmission        - Log10 of mission Pc                                
%
%   L10PremSearch       - Array of Log10 of Prem values               [1xN]
%
%   L10Pcummd           - Array of Log10 of median Pcum values        [1xN]
%
%   L10Pcumlo           - Array of Log10 of -1 sigma Pcum values      [1xN]
%
%   L10Pcumhi           - Array of Log10 of +1 sigma Pcum values      [1xN]
% 
%   Nmvrmd              - Array of median number of maneuvers         [1xN]
%
%   Nmvrlow             - Array of -1 sigma number of maneuvers       [1xN]
%
%   Nmvrhi              - Array of +1 sigma number of maneuvers       [1xN]
%
% =========================================================================
%
% Output:
%
%   out                 - structure containing interpolated median  
%                         and +/- 1-sigma values for Log10 of Prem,  
%                         achievednumber of maneuvers, and achieved 
%                         Pcum
%
% =========================================================================
%
% Dependencies:
%
%   None
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

out.L10Premmd = interp_LogPrem( ...
    L10Pcummd,L10PremSearch,L10Pcmission);
out.L10Premlo = interp_LogPrem( ...
    L10Pcumhi,L10PremSearch,L10Pcmission);
out.L10Premhi = interp_LogPrem( ...
    L10Pcumlo,L10PremSearch,L10Pcmission);

out.NmvrmdAchieved = interp1(L10PremSearch,Nmvrmd,out.L10Premmd);
out.NmvrloAchieved = interp1(L10PremSearch,Nmvrlo,out.L10Premmd);
out.NmvrhiAchieved = interp1(L10PremSearch,Nmvrhi,out.L10Premmd);

out.L10PcummdAchieved = interp1(L10PremSearch,L10Pcummd,out.L10Premmd);
out.L10PcumloAchieved = interp1(L10PremSearch,L10Pcumlo,out.L10Premlo);
out.L10PcumhiAchieved = interp1(L10PremSearch,L10Pcumhi,out.L10Premhi);

temp = [out.L10PcumloAchieved out.L10PcumhiAchieved];
out.L10PcumloAchieved = min(temp);
out.L10PcumhiAchieved = max(temp);

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
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
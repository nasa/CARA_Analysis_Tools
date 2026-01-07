function LogPremint = interp_LogPrem(LogPcum,LogPrem,LogPcumgoal)
% interp_LogPrem - Interpolate LogPrem values
%
% Syntax: LogPremint = interp_LogPrem(LogPcum,LogPrem,LogPcumgoal)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Interpolate LogPrem values. Returns a Prem value
%              corresponding to the desired Pcum through linear 
%              interpolation of Prem and corresponding Pcum values
%
% =========================================================================
%
% Input:
%
%   LogPcum      - Log of Pcum values                                 [1xN]
%
%   LogPrem      - Log of Prem values corresponding to LogPcum values [1xN]
%
%   LogPcumgoal  - Desired LogPcum value for which to interpolate a
%                  corresponding LogPrem
% =========================================================================
%
% Output:
%
%   LogPremint   - Interpolated Prem value corresponding to LogPcumgoal
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

if LogPcum(1) > LogPcumgoal
    
    LogPremint = LogPrem(1);
    
elseif LogPcum(end) < LogPcumgoal
    
    LogPremint = min(LogPcumgoal,LogPrem(end));
    
else

    n1 = max(find((LogPcum <= LogPcumgoal)));
    n2 = n1+1;
    LogPremint = interp1(LogPcum(n1:n2),LogPrem(n1:n2), ...
                         LogPcumgoal,'linear','extrap');
                     
    LogPremint = min(LogPcumgoal,LogPremint);
    
end

return;
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
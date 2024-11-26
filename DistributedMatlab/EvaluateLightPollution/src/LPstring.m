function LP = LPstring(Impact,LowToHigh)
% LPstring - Convert light pollution (LP) impact level to an LP impact string
% Syntax: LP = LPstring(Impact,LowToHigh);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

if isnan(Impact) || Impact < 0
    LP = 'N/A';
elseif Impact == 0
    LP = 'None';
elseif Impact < LowToHigh(1)
    LP = 'Very Low';
elseif Impact < LowToHigh(2)
    LP = 'Low';
elseif Impact < LowToHigh(3)
    LP = 'Medium';
elseif Impact < LowToHigh(4)
    LP = 'High';
else
    LP = 'Very High';
end

LP = upper(LP);

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall    |  2024-Nov  | Initial version.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
function tag = make_timetag(timestring)
% make_timetag - Alters a timestring into a format used by EventRate
%                functions
%
% Syntax: tag = make_timetag(timestring)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Alters a timestring into a format used by EventRate
%              functions by removing '-' and ':' characters and truncating
%              the year
%
% =========================================================================
%
% Input:
%
%   timestring   - Date string formatted as:
%
%                   'YYYY-MM-DD HH-mm-SS'
%
% =========================================================================
%
% Output:
%
%   tag          - Date string formatted as:
%
%                   'YYMMDD_HHmmSS'
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
tag = timestring(3:19);
tag = strrep(tag,'-','');
tag = strrep(tag,':','');
tag = strrep(tag,' ','_');
tag = strrep(tag,'T','_');

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
function tag = make_timetag(timestring)
% make_timetag - Reformat timestring into timetag.
% Syntax: make_timetag(timestring);
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

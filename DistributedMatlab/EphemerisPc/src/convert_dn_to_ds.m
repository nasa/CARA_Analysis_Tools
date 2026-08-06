function [ds,dscompact] = convert_dn_to_ds(dn)
% convert_dn_to_ds - Convert a Matlab date number to a date string,
%                    trimming to the nearest second if appropriate, and
%                    also making a compact version appropriate for using
%                    in file names.
%
% Syntax: [ds, dscompact] = convert_dn_to_ds(dn);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    dn         -   Matlab date number
%
% =========================================================================
%
% Output:
%
%   ds          -   Date string (yyyy-mm-dd HH:MM:SS.FFF,
%                   or without .FFF if fractional seconds
%                   are zero)
%
%   dscompact   -   Compact date string (yyyymmdd_HHMMSS)
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Full date string
ds = datestr(dn,'yyyy-mm-dd HH:MM:SS.FFF');

% Date string trimmed to eliminate all-zero .FFF fields
if strcmpi(ds(20:23),'.000')
    ds = ds(1:19);
end

% Date string made compact into yyyymmdd_HHMMSS format, rounding .FFF down
dscompact = ds(1:19);
dscompact = strrep(dscompact,':','');
dscompact = strrep(dscompact,'-','');
dscompact = strrep(dscompact,' ','_');

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
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
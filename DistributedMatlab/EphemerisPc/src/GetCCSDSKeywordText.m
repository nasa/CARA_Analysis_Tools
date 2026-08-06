function str = GetCCSDSKeywordText(data,keyword)
% GetCCSDSKeywordText - Get a keyword string from CCSDS OEM file text
%                       data.
%
% Syntax: str = GetCCSDSKeywordText(data, keyword);
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
%    data       -   CCSDS OEM file text data cell array   {Nx1}
%
%    keyword    -   Keyword to search for
%
% =========================================================================
%
% Output:
%
%   str         -   Keyword value string
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

ndx = find(~isempty_cell(regexpi(data,keyword)),1,'first');
str = textscan(char(data(ndx)),'%*s %*s %s');
str = str{1}; str = str{1};

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
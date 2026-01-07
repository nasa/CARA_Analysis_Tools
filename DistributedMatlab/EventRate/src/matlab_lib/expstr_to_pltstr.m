function pltstr = expstr_to_pltstr(expstr)
% expstr_to_pltstr - Convert a numerical string to Matlab plot notation
%
% Syntax: pltstr = expstr_to_pltstr(expstr)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Convert a numerical string to Matlab plot notation.
%              Targets strings like '1e-3' and '2.5e-5' and converts  
%              them to '10^{-3}' and '2.5{\times}10^{-5}'
%
% =========================================================================
%
% Input:
%
%   expstr   - input string with numerical notation              
% =========================================================================
%
% Output:
%
%   phat     - reformatted string with plot notation
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
expstr = lower(strtrim(expstr));

k = strfind(expstr,'e');

if ~isempty(k) && (k > 1) && (k < length(expstr))
    
    % Modify exp notation to Matlab plot notation
    
    if (k == 2) && strcmpi(expstr(1),'1')
        % Change 1e-3 into 10^{-3}
        pltstr = ['10^{' expstr(3:end) '}'];
    else
        % Change 2.5e-5 into 2.5{\times}10^{-5}
        pltstr = [expstr(1:k-1) '{\times}10^{' expstr(k+1:end) '}'];
    end
    
else
    
    % Use original string as plot string
    
    pltstr = expstr;
    
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
function M = SatCon1Limit(hkm)
% SatCon1Limit - Return brightness limit using formula from SatCon-1 Workshop.
% Syntax: M = SatCon1Limit(hkm);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%    hkm - height (km)
%
% =========================================================================
%
% Output:
%    M   - Magnitude limit  
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% SatCon-1 Workshop Report brightness limit
M = 7+2.5*log10(hkm/550);
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
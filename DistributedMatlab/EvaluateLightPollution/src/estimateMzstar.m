function Mzstar = estimateMzstar(Fstar,Mmd,Flo,Mlo)
% estimateMzstar - Calculate Mzstar
% Syntax: Mzstar = estimateMzstar(Fstar,Mmd,Flo,Mlo)
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

% Calculate the Gaussian sigma value either above or below the median
sqrt2 = sqrt(2);
sig = (Mlo-Mmd)/sqrt2/erfinv(-1+2*Flo);

Mzstar = Mmd-sqrt2*sig*erfinv(1-2*Fstar);

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
function [phat,pci] = ConjDist_binofit(x,n)
% ConjDist_binofit - Call binofit, but catch the special case of n=0.
%
% Syntax: [phat,pci] = ConjDist_binofit(x,n)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Call binofit, but catch the special case of n=0.
%
% =========================================================================
%
% Input:
%
%   x      - Number of successes                                      [Nx1]
%
%   n      - Number of trials. The i-th entry of n must be greater    [Nx1]
%            than or equal to the i-th entry of x                     
% =========================================================================
%
% Output:
%
%   phat   - MLE of the progability of success for a binomial         [Nx1] 
%            distribution with x successes over n trials
%
%   pci    - 95% confidence interval lower/upper bounds               [Nx2] 
%            corresponding to each estimate phat
%
% =========================================================================
%
% Dependencies:
%
%   None
%
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
if (n > 0)
    [phat,pci] = binofit(x,n);
elseif (n == 0)
    phat = 0; pci = [0 0];
else
    phat = NaN; pci = [NaN NaN];
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
function [amd,alo,ahi] = ConjDist_set_resample(nset,xmd,xlo,xhi, ...
                                               amd,alo,ahi, ...
                                               Nsamp,alpha,limits)
% ConjDist_resample - Calculate the average for a set of numbers, and 
% resample to get the 1-alpha confidence range. Checks for similarity with
% full data set.
%
% Syntax: [amd,alo,ahi] = ConjDist_set_resample(nset,xmd,xlo,xhi, ...
%                                               amd,alo,ahi, ...
%                                               Nsamp,alpha,limits)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculate the average for a set of numbers, and resample to 
% get the 1-alpha confidence range. Also, if this is the SVI (nset=2) or 
% nonSVI (nset=3) set, then check if it is identical to the ALL (nset=1) 
% set, and if so, use those results.
%
% =========================================================================
%
% Input:
%
%   nset   - set index to process (equal to 1, 2, or 3)
%
%   xmd    - Median values                                            [3xN]
%
%   xlo    - low-deviation (-1 sigma) values                          [3xN]
%
%   xhi    - high-deviation (+1 sigma) values                         [3xN]
%
%   amd    - average median value for each set                        [3x1]
%
%   alo    - avg low-deviation (-1 sigma) conf value for each set     [3x1]
%
%   ahi    - avg high-deviation (+1 sigma) conf value for each set    [3x1]
%
%   Nsamp  - Optional - Number of points to resample 
%            (default = 40)
%
%   alpha  - Optional - Confidence interval range 
%            (default = 0.05 for 95% confidence interval)
%
%   limits - Optional - Lower/upper limits for resampled data         [1x2]
%            (default = [-inf, inf] (no limits))
% =========================================================================
%
% Output:
%
%   amd    - average median value for each set                        [3x1]
%
%   alo    - avg -1 sigma confidence value for each set               [3x1]
%
%   ahi    - avg +1 sigma confidence value for each set               [3x1]
%
% =========================================================================
%
% Dependencies:
%
%   ConjDist_resample.m
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

% Initializations

Nargin = nargin;
if (Nargin < 10); limits = [-Inf,Inf]; end
if (Nargin <  9); alpha = 0.05; end
if (Nargin <  8); Nsamp = 40; end

% If this is not the first set, check if it is the same as the first set

if (nset > 1)                     && ...
    isequal(xmd(nset,:),xmd(1,:)) && ...
    isequal(xlo(nset,:),xlo(1,:)) && ...
    isequal(xhi(nset,:),xhi(1,:))
    % This set same as first set, so use those results
    amd(nset) = amd(1);
    alo(nset) = alo(1);
    ahi(nset) = ahi(1);
else
    % Resample this set
    [amd(nset),alo(nset),ahi(nset)] = ...
        ConjDist_resample(xmd(nset,:), ...
                          xlo(nset,:), ...
                          xhi(nset,:), ...
                          Nsamp,alpha,limits);
end

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
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
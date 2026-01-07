function [xmed] = median_sorted(xsort)
% median_sorted - Median from sorted data
%
% Syntax: [xmed] = median_sorted(xsort)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Median from sorted data
%
% =========================================================================
%
% Input:
%
%   xsort - Sorted data                                               [1xN]
%
% =========================================================================
%
% Output:
%
%   xmed  - median
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

% Number of elements in sorted data set

Nx = numel(xsort);

% Median index

nmd = 1+(Nx-1)/2;

% Handle even/odd element numbers appropriately

if mod(Nx,2) == 0
    % Median for even number of samples
    xmed = 0.5*(xsort(floor(nmd))+xsort(ceil(nmd)));
else
    % Median for odd number of samples
    xmed = xsort(nmd);
end

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
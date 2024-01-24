function [xmed] = median_sorted(xsort)

% Median from the sorted data

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
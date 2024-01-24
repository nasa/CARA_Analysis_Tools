function [xlo,xhi] = conf_range_sorted(xsort,alpha)

% Confidence range from sorted data

% Number of elements in sorted data set

Nx = numel(xsort);

% Half alpha value

h = alpha/2;

% Low and high limits

nlo = min(Nx,max(1,floor(Nx*h)));
xlo = xsort(nlo);

nhi = min(Nx,max(1,ceil(Nx*(1-h))));
xhi = xsort(nhi);

return;
end
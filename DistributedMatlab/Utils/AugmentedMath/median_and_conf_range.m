function [xmd,xlo,xhi] = median_and_conf_range(x,alpha)

% Calculate the median value and confidence range for a distribution

% Initializations and defaults

Nargin = nargin;
if Nargin < 2; alpha = []; end
if isempty(alpha); alpha = 0.05; end

% Calculate median and range

if isempty(x)
    % Return null values for null input
    xmd = []; xlo = []; xhi = [];
else
    % Calculate median and range
    xsort = sort(x);
    [xmd] = median_sorted(xsort);
    [xlo,xhi] = conf_range_sorted(xsort,alpha);
end

return
end
function [amd,alo,ahi] = ConjDist_resample(xmd,xlo,xhi,Nsamp,alpha,limits)
% ConjDist_resample - Calculate the average for a set of numbers, and 
% resample to get the 1-alpha confidence range
%
% Syntax: [amd,alo,ahi] = ConjDist_resample(xmd,xlo,xhi,Nsamp,alpha,limits)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculate the average for a set of numbers, and 
% resample to get the 1-alpha confidence range
%
% =========================================================================
%
% Input:
%
%   xmd    - Median values                                            [Nx1]
%
%   xlo    - low-deviation (-1 sigma) values                          [Nx1]
%
%   xhi    - high-deviation (+1 sigma) values                         [Nx1]
%
%   Nsamp  - Optional - number of points to resample
%            (default = 40)
%
%   alpha  - Optional - confidence interval range
%            (default = 0.05: 1 - alpha = 95% confidence interval)
%
%   limits - Optional - Lower/upper limits for resampled data         [1x2]
%            (default = [-inf, inf])
% =========================================================================
%
% Output:
%
%   amd    - average median value
%
%   alo    - average low-deviation (-1 sigma) confidence value
%
%   ahi    - average high-deviation (+1 sigma) confidence value
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

% Initializations

Nargin = nargin;
if (Nargin < 6); limits = [-Inf,Inf]; end
if (Nargin < 5); alpha = 0.05; end
if (Nargin < 4); Nsamp = 40; end

% Median

amd = mean(xmd);

% Calculate 1-sigma low and high deviations

chifact = sqrt(chi2inv(1-alpha,1));
siglo = (xmd-xlo)/chifact;
sighi = (xhi-xmd)/chifact;

% Calculate multiple samples per x point

Nx = numel(xmd);
Nxs = Nx*Nsamp;

xsamp = zeros(Nx,Nsamp);

for nx=1:Nx
    rsamp = randn(1,Nsamp);
    ndx = (rsamp < 0); % Low deviations
    xsamp(nx,ndx) = xmd(nx)+siglo(nx)*rsamp(ndx);
    ndx = ~ndx; % High deviations
    xsamp(nx,ndx) = xmd(nx)+sighi(nx)*rsamp(ndx);
end        

xsrt = sort(xsamp(:));

xsrt(xsrt < limits(1)) = limits(1);
xsrt(xsrt > limits(2)) = limits(2);

al = alpha/2;
n1 = max(floor(al*Nxs),1);
n2 = min(ceil(Nxs-al*Nxs),Nxs);

alo = xsrt(n1);
ahi = xsrt(n2);

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
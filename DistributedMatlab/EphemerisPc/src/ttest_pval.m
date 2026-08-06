function [p,t,dfe] = ttest_pval(xb,sx,nx,yb,sy,ny)
% ttest_pval - T-test p value (two-tailed test for unequal variances).
%
% Syntax: [p, t, dfe] = ttest_pval(xb, sx, nx, yb, sy, ny);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
%    xb         -   x mean
%
%    sx         -   x standard deviation
%
%    nx         -   x size
%
%    yb         -   y mean
%
%    sy         -   y standard deviation
%
%    ny         -   y size
%
% =========================================================================
%
% Output:
%
%   p           -   Two-tailed p value
%
%   t           -   T statistic
%
%   dfe         -   Degrees of freedom (Welch-Satterthwaite)
%
% =========================================================================
%
% References:
%
%   See Matlab ttest2.m
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

difference = xb-yb;
s2x = sx.^2; s2xbar = s2x ./ nx;
s2y = sy.^2; s2ybar = s2y ./ ny;
se = sqrt(s2xbar+s2ybar);
t = difference/se;
dfe = (s2xbar + s2ybar).^2 ./ (s2xbar.^2 ./ (nx-1) + s2ybar.^2 ./ (ny-1));
p = 2 * tcdf(-abs(t),dfe);

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
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
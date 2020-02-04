function [g,converged] = k2b_cont_frac(q,b_tolerance,max_iterations)
% =========================================================================
%
% Shepperd (1985) continued fraction iteration.
%
% =========================================================================
%
% REFERENCES:
%
%    S. W. Shepperd, "Universal Keplerian State Transition Matrix"
%    Celestial Mechanics, 35, 129-144, 1985.  [Hereafter S85]
%
% =========================================================================
%
%
% =========================================================================

% See S85 equations (A.20)-(A.33).

n = 0;
l = 3;
d = 15;
k = -9;
a = 1;
b = 1;
g = 1;

converged = true;
niter = 0;

while (abs(b) > b_tolerance)
    
    niter = niter+1;

    k = -k;
    l = l + 2;
    d = d + 4 * l;
    n = n + (1 + k) * l;
    a = d / (d - n * a * q);
    b = (a - 1) * b;
    g = g + b;
    
    % disp(['n = ' num2str(niter) '  b = ' num2str(b)]);
    
    if (niter > max_iterations)
        converged = false;
        break;
    end
    
end

% disp(['q,g,i = ' num2str(q,'%0.18e') ' ' num2str(g,'%0.18e') ' ' num2str(niter)]);

return;
end


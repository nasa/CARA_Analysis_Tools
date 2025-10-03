function [c,tinc,Finc,rankAmat,x] = TimeParabolaFit(t,F)
% TimeParabolaFit - Calculate best-fit parabola to a function using only
%                   (t,F) points nearest the minimum F point.
%
% Syntax: [y,tinc,Finc,rankAmat,x] = TimeParabolaFit(t,F);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% DESCRIPTION:
%
%   Selecting the fewest points possible from the t and F arrays passed in,
%   calculate coefficients for a best-fit parabola using the minimum F
%   point as a basis.
%
%   The "c" output contains coefficients of the parabola such that an array
%   of x-values spanning "t" would fit the "F" values used as best as
%   possible using the following equation:
%     y = c(1).*x.^2 + c(2).*x + c(3);
%
% =========================================================================
%
% INPUT:
%
%   t = Array of times (or x-values) associated with the parabola. [1xn]
%   F = Array of function values at "t" times. Must be the same size as the
%       "t" array. [1xn]
%
% =========================================================================
%
% OUTPUT:
%
%   c = Coefficients of the best-fit parabola [3x1]
%   tinc = Times from the "t" array used for the best-fit parabola solution
%   Finc = Values from the "F" array used for the best-fit parabola
%          solution
%   rankAmat = Rank of the design matrix used in creating the solution
%   x = Coefficients of the best-fit parabola using normalized times [3x1]
%
% =========================================================================
%
% Initial version: Mar 2023; Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

% Restrict to unique times, to avoid matrix singularities
[t,indx] = unique(t); F = F(indx);

% Number of time points
Nt = numel(t);
if Nt < 3
    error('Minimum of three unique points required for parabolic fit');
end

% Sort all points by ascending F value
[Fvec,srt] = sort(F);
tvec = t(srt);

% Make column vectors
tvec = reshape(tvec,[Nt 1]);
Fvec = reshape(Fvec,[Nt 1]);

% Initially include the 3 unique points with smallest F values
inc = false(size(tvec)); inc(1:3) = true;
tinc = tvec(inc);

% Solve using centered, normalized time which improves matrix conditioning
tmin = min(tinc);
tmax = max(tinc);
tdel = 0.5*(tmax-tmin);
tmid = 0.5*(tmax+tmin);
z = (tinc-tmid)/tdel;

% Initial design matrix for 3-point solution, and associated rank
Amat = [z.^2 z ones(size(z))];
rankAtol = 1000*max(size(Amat))*eps(norm(Amat));
rankAmat = rank(Amat,rankAtol);

% Add more points if required to ensure full rank matrix
while rankAmat < 3 && any(~inc)
    fnd = find(~inc,max(1,3-sum(inc)));
    inc(fnd) = true;
    tinc = tvec(inc);
    tmin = min(tinc);
    tmax = max(tinc);
    tdel = 0.5*(tmax-tmin);
    tmid = 0.5*(tmax+tmin);
    z = (tinc-tmid)/tdel;    
    Amat = [z.^2 z ones(size(z))];
    rankAtol = 1000*max(size(Amat))*eps(norm(Amat));    
    rankAmat = rank(Amat,rankAtol);
end

% Psuedoinverse solution parabola coefficients as a function of z
Finc = Fvec(inc);
x = pinv(Amat)*Finc;

% Convert to coefficients of t instead of z
trat = tmid/tdel;
tratx1 = trat*x(1);
c = [x(1)/tdel^2; ...
    (x(2)-2*tratx1)/tdel; ...
    trat*(tratx1-x(2))+x(3)];

return;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2023-Mar-01 | Initial Development
% L. Baars       | 2025-Aug-06 | Updated documentation for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

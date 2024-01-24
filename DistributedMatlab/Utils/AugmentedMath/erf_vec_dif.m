function d = erf_vec_dif(a,b)
%
% Calculate the difference d = erf(a) - erf(b), but use erfc for cases of
% large positive or negative values of a and b, which provides improved
% accuracy.
%

% % Check input
% Sa = size(a); Na = prod(Sa);
% if ~isequal(Sa,size(b))
%     error('Invalid input - unequal (a,b) dimensions');
% end
% 
% % Make both column vectors
% a = reshape(a,[Na 1]);
% b = reshape(b,[Na 1]);

% Min and max
% ab = [a; b];
ab = [a(:) b(:)];
minab = min(ab,[],2);
maxab = max(ab,[],2);

% Large measure for argument of erf
large = 3;

% Initialize array of difference d = erf(a) - erf(b)
% d = NaN([Na 1]);
d = NaN(size(a));

% For large positive a & b values, use erfc for better accuracy
set1 = minab >  large;
d(set1) = erfc(b(set1)) - erfc(a(set1));

% For large negative a & b values, use erfc with negated arguments
set2 = maxab < -large;
d(set2) = erfc(-a(set2)) - erfc(-b(set2));

% Use erf for all other cases
set3 = ~set1 & ~set2;
d(set3) = erf(a(set3)) - erf(b(set3));

% % Reshape for output
% d = reshape(d,Sa);

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% D.Hall         | 2020-12-17 |  Implemented to improve the accuracy
%                                of the function Pc2D_SquareErrFunc.m for
%                                very small Pc values, enabling it to
%                                efficiently reproduce the results from
%                                the numerical integration function
%                                Pc2D_RotatedSquare.m.
% D.Hall         | 2021-12-23 |  Vectorized algorithm.
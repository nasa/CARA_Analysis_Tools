function d = erf_dif(a,b)
%
% Calculate the difference d = erf(a)-erf(b), but use erfc for cases of
% large positive or negative values of a and b, which provides improved
% accuracy.
%

% Large measure for argument of erf
large = 3;

% Handle different (a,b) cases
if min(a,b) > large
    % For large positive a & b values, use erfc for better accuracy
    d = erfc(b)-erfc(a);
elseif max(a,b) < -large
    % For large negative a & b values, use erfc with negated arguments
	d = erfc(-a)-erfc(-b);
else
    % Use erf for all other cases
    d = erf(a)-erf(b);
end

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
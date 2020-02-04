function Csym = cov_make_symmetric(C)
%
% cov_make_symmetric - Make a covariance matrix diagonally symmetric if
%                      required.
%
% Syntax: Csym = cov_make_symmetric(C)
%
% Inputs:
%   C    - Input covariance matrix, must be [NxN]
%
% Outputs:
%   Csym - Symmetrized version of covariance matrix [NxN]
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
% See also: None
%
% August 2019; Last revision: 2019-AUG-28
%
% ----------------- BEGIN CODE -----------------

% Check for bad covariance matrix input
szC = size(C);
if (numel(szC) ~= 2)
    error('Array needs to be a 2D matrix to make symmetric.');
end
if (szC(1) ~= szC(2))
    error('Matrix needs to be square to make symmetric.');
end

% Calculate transpose
Ct = C';

% Check existing status of diagonal symmetry
if isequal(C,Ct)
    
    % Original matrix is already diagonally symmetric
    Csym = C;
    
else
    
    % Average out any off-diagonal asymmetries
    Csym = (C+Ct)/2;
    
    % Reflect about diagonal to ensure diagonal symmetry absolutely
    Csym = triu(Csym,0)+triu(Csym,1)';
    
end

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
% D.Hall         | 2019-AUG-28 | Initial Development
%

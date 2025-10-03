function [TF] = isclose(a, b, tolLevel)
% isclose - This function determines if two arrays are element-wise equal
%           within a tolerance.
%
% Syntax: TF = isclose(a, b, tolLevel);
%         TF = isclose(a, b);
%
% =========================================================================
%
% Description:
%
%   The algorithm is adapted from the "Relative and absolute tolerances"
%   section of the Matlab isapprox documentation. The tolerance values are
%   based on the tolerance levels for non-single data types within the same
%   document.
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    a - Input array of doubles
%
%    b - Input array of doubles, must be the same size as the "a" variable
%
%    tolLevel - (Optional) Tolerance value of the comparison between "a"
%               and "b". Valid values are:
%                 'verytight' - tolerance of 1e-15
%                 'tight'     - tolerance of 1e-12
%                 'loose'     - tolerance of 1e-8
%                 'veryloose' - tolerance of 1e-3
%               Defaults to 'verytight'
%
% =========================================================================
%
% Output:
%
%    TF - Boolean array the same size as "a" and "b". True values indicate
%         that the corresponding values in the "a" and "b" arrays are equal
%         to within the specified tolerance.
%
% =========================================================================
%
% Initial version: Sep 2025;  Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------

    if nargin == 2
        tolLevel = 'verytight';
    elseif nargin ~= 3
        error('Incorrect number of parameters passed in');
    end

    if strcmpi(tolLevel,'verytight')
        absTol = 1e-15;
        relTol = 1e-15;
    elseif strcmpi(tolLevel, 'tight')
        absTol = 1e-12;
        relTol = 1e-12;
    elseif strcmpi(tolLevel, 'loose')
        absTol = 1e-8;
        relTol = 1e-8;
    elseif strcmpi(tolLevel, 'veryloose')
        absTol = 1e-3;
        relTol = 1e-3;
    else
        error('Incorrect tolLevel passed in, valid values are ''verytight'', ''tight'', ''loose'', or ''veryloose''');
    end

    TF = abs(a - b) <= max(absTol, relTol .* max(abs(a),abs(b)));
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% L. Baars       | 09-16-2025 | Initial Development

% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

function [Constant,VarVal] = isconstant(x,VarTol)
% isconstant - Evaluate if all of the elements of an array x are constant
% to within a specified tolerance
%
% Syntax: [Constant,VarVal] = isconstant(x,VarTol)
%
% Inputs:
%   x         - Array to be analyzed (any dimension)
%   VarTol    - Tolerance to evaluate x element variations (default: 1e-4)
%
% Outputs:
%   Constant - Logical indicating if all x array elements are constant to
%              within the specified variation tolerance level VarTol:
%                 true  => x values constant to within tolerance
%                 false => x values vary beyond tolerance
%   VarVal   - Variation value used to evaluate constance
%
% Example/Validation Cases:
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: none
%
% August 2019; Last revision: 2019-SEP-04
%
% ----------------- BEGIN CODE -----------------

    % Set default variation tolerance
    if (nargin < 2) || isempty(VarTol)
        VarTol = 1e-4;
    end
    
    % Evaluate x array for variations
    Nx = numel(x);
    if Nx <= 1
        % Arrays with one or zero elements can have no variations,
        % so are evaluated as constant
        Constant = true;
        VarVal = 0;
    else
        % Min and max of linear-indexed x array
        minx = min(x(:));
        maxx = max(x(:));
        % Evaluate variation in x array elements
        if minx == maxx
            % Equal min and max values indicate no variations
            Constant = true;
            VarVal = 0;
        else
            % Compare difference between unequal min and max x values to 
            % cutoff value defined by specified tolerance
            difx = maxx-minx;
            medx = 0.5*(minx+maxx);
            valx = max([difx abs(medx) abs(minx) abs(maxx)]);
            Constant = (difx <= valx*VarTol);
            VarVal = difx/valx;
        end
    end
    
    return
    
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D.Hall         | 2019-SEP-04 | Initial Development
%

function A = AspectAvgProjArea(H,W,L)
% A - Aspect-averaged projected area of a box = surface-area/4
% Syntax: A = AspectAvgProjArea(H,W,L)
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    H          -   Height
%    W          -   Width
%    L          -   Length
%
% =========================================================================
%
% Output:
%
%    a          -   Aspect-averaged projected area of a box   
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------


% Aspect-averaged projected area of a box = surface-area/4
A = 0.5*(H*W+H*L+W*L);
return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall    |  2024-Nov  | Initial version.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
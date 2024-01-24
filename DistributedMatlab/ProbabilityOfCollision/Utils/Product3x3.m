function [Out] = Product3x3(a,b)
% Product3x3 - vectorized 3x3 matrix multiplication routine
%
% Syntax: [Out] = Product3x3(a,b);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   a       -   Left matrix                                           [nx9]
%   b       -   Right matrix                                          [nx9]
%
% =========================================================================
%
% Output:
%
%   Out     -   Matrix product                                        [nx9]
%
% =========================================================================
%
% Notes:
%
%   Inputs and outputs are formatted as [col1, col2, col3]:
%   [(1,1) (2,1) (3,1) (1,2) (2,2) (3,2) (1,3) (2,3) (3,3)]
%
% =========================================================================
%
% Initial version: May 2022;  Latest update: Jul 2023
%
% ----------------- BEGIN CODE -----------------

    Out = [a(:,1).*b(:,1)+a(:,2).*b(:,4)+a(:,3).*b(:,7) ...
           a(:,1).*b(:,2)+a(:,2).*b(:,5)+a(:,3).*b(:,8) ...
           a(:,1).*b(:,3)+a(:,2).*b(:,6)+a(:,3).*b(:,9) ...
           a(:,4).*b(:,1)+a(:,5).*b(:,4)+a(:,6).*b(:,7) ...
           a(:,4).*b(:,2)+a(:,5).*b(:,5)+a(:,6).*b(:,8) ...
           a(:,4).*b(:,3)+a(:,5).*b(:,6)+a(:,6).*b(:,9) ...
           a(:,7).*b(:,1)+a(:,8).*b(:,4)+a(:,9).*b(:,7) ...
           a(:,7).*b(:,2)+a(:,8).*b(:,5)+a(:,9).*b(:,8) ...
           a(:,7).*b(:,3)+a(:,8).*b(:,6)+a(:,9).*b(:,9)];
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% L. Baars       | 05-19-2022 | Initial development
% E. White       | 07-12-2023 | Added compliant documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

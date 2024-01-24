classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Product3x3_UnitTest < matlab.unittest.TestCase
% Product3x3_UnitTest - Unit test for Product3x3
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Jul 2023;  Latest update: Jul 2023
%
% ----------------- BEGIN CODE -----------------
    
    methods (Test)
        % Tests 1000 random 3x3 matrices
        function testDiagonal (testCase)
            A = 2000 * rand(1000, 9) - 1000;
            B = 2000 * rand(1000, 9) - 1000;
            
            P = Product3x3(A, B);
            
            for i = 1:1000
                M = reshape(P(i, :), 3, 3)';
                T = reshape(A(i, :), 3, 3)' * reshape(B(i, :), 3, 3)';
                
                testCase.verifyEqual(M, T);
            end
        end
    end
    
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% E. White       | 07-03-2023 | Initial development

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

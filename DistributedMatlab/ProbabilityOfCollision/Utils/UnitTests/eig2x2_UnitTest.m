classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        eig2x2_UnitTest < matlab.unittest.TestCase
% eig2x2_UnitTest - Unit test for eig2x2
%
% =========================================================================
%
% Copyright (c) 2023-2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Jul 2023;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------
    
    methods (Test)
        % Tests 1000 random 2x2 real diagonal matrices
        function testDiagonal (testCase)
            a = 2000 * rand(1000, 1) - 1000;
            b = zeros(1000, 1);
            d = 2000 * rand(1000, 1) - 1000;
            
            [V1, V2, L1, L2] = eig2x2([a, b, d]);
            
            for i = 1:1000
                [V, D] = eig([a(i), b(i); b(i), d(i)], 'vector');
                [D, idx] = sort(D, 'descend');
                V = V(:, idx);
                testCase.verifyEqual([L1(i); L2(i)], D, 'RelTol', 1E-12);
                testCase.verifyEqual([V1(i, :)', V2(i, :)'], V);
            end
        end
        
        % Tests 2x2 zero matrix
        function testZeroMatrix (testCase)
            [V1, V2, L1, L2] = eig2x2([0, 0, 0]);
            
            testCase.verifyEqual([L1 L2], zeros(1, 2));
            testCase.verifyEqual([V1; V2], eye(2));
        end
        
        % Tests 1000 random 2x2 real symmetric matrices against MATLAB's
        % eig() function
        
        % Note that the eigenvectors sometimes differ by a factor of -1
        function testRandomMatrices (testCase)
            a = 2000 * rand(1000, 1) - 1000;
            b = 2000 * rand(1000, 1) - 1000;
            d = 2000 * rand(1000, 1) - 1000;
            
            [V1, V2, L1, L2] = eig2x2([a, b, d]);
            
            for i = 1:1000
                [V, D] = eig([a(i), b(i); b(i), d(i)], 'vector');
                [D, idx] = sort(D, 'descend');
                V = V(:, idx);
                testCase.verifyEqual([L1(i); L2(i)], D, 'RelTol', 1E-5, 'AbsTol', 1E-5);
                f = abs([V1(i, :)', V2(i, :)'] ./ V) - 1;
                f(f <= 1E-5) = 0;
                testCase.verifyEqual(f, zeros(2));
            end
        end
        
        % Tests 1000 random 2x2 real symmetric matrices with small
        % off-diagonal values against MATLAB's eig() function
        
        % Note that the eigenvectors sometimes differ by a factor of -1
        function testSmallOffDiagonalMatrices (testCase)
            a = 2000 * rand(1000, 1) - 1000;
            b = 2E-2 * rand(1000, 1) - 1E-2;
            d = 2000 * rand(1000, 1) - 1000;
            
            [V1, V2, L1, L2] = eig2x2([a, b, d]);
            
            for i = 1:1000
                [V, D] = eig([a(i), b(i); b(i), d(i)], 'vector');
                [D, idx] = sort(D, 'descend');
                V = V(:, idx);
                testCase.verifyEqual([L1(i); L2(i)], D);
                testCase.verifyEqual([V1(i, :)', V2(i, :)'], V);
            end
        end
        
        % Test the default covariance, which can cause the determinant
        % calculation to fail because the diagonals are so much larger than
        % the off-diagonals that floating point subtraction of the numbers
        % fails
        function testDefaultCovariance (testCase)
            a = 4.06806226869326435234562435e+15;
            b = 251651.25;
            d = 4.06806226869326435234562435e+15;
            
            [V1, V2, L1, L2] = eig2x2([a, b, d]);
            
            [V, D] = eig([a, b; b, d],'vector');
            [D, idx] = sort(D, 'descend');
            V = V(:, idx);
            testCase.verifyEqual([L1; L2], D);
            testCase.verifyEqual([V1(1, :)', V2(1, :)'], V);
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
% E. White       | 07-12-2023 | Initial development
% L. Baars       | 11-29-2024 | Updated testDiagonal unit test to check
%                               against eig() outputs. Added unit test for
%                               the default covariance.

% =========================================================================
%
% Copyright (c) 2023-2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        eig2x2_UnitTest < matlab.unittest.TestCase
% eig2x2_UnitTest - Unit test for eig2x2
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
        % Tests 1000 random 2x2 real diagonal matrices
        function testDiagonal (testCase)
            a = 2000 * rand(1000, 1) - 1000;
            b = zeros(1000, 1);
            d = 2000 * rand(1000, 1) - 1000;
            
            [V1, V2, L1, L2] = eig2x2([a, b, d]);
            
            % Swap values of a and d so a > d (for checking eigenvalues)
            t = a;
            i = a < d;
            a(i) = d(i);
            d(i) = t(i);
            
            testCase.verifyEqual(a, L1, 'RelTol', 1E-5, 'AbsTol', 1E-5);
            testCase.verifyEqual(d, L2, 'RelTol', 1E-5, 'AbsTol', 1E-5);
            testCase.verifyEqual(V1, repmat([1, 0], 1000, 1));
            testCase.verifyEqual(V2, repmat([0, 1], 1000, 1));
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

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../../AugmentedMath')}) ...
        ECI2RIC_UnitTest < matlab.unittest.TestCase
% ECI2RIC - Unit test for ECI2RIC
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Aug 2023;  Latest update: Oct 2025
%
% ----------------- BEGIN CODE -----------------
    
    properties (TestParameter)
        
        relTol = {1e-12};
        
        defaultR = {[378.39559 4305.721887 5752.767554]};
        defaultV = {[2.360800244 5.580331936 -4.322349039]};

        three = {[ 44.5757544811362   81.6751751052616 -67.8687662707124;
                   81.6751751052616  158.453402956163  -128.616921644857;
                  -67.8687662707124 -128.616921644858   105.490542562701]};
        threeOut = {[ 0.6597  -5.963   1.253;
                     -5.963    305    -19.41;
                      1.253   -19.41   2.86]};
              
        six = {[ 3.748e-11   2.066e-11  -9.1357e-11 -3.0881e-12 -4.2085e-13  1.2413e-11; ...
                 2.066e-11   2.4238e-11 -4.8871e-11 -2.3602e-12  2.1508e-13 -4.9476e-12; ...
                -9.1357e-11 -4.8871e-11  3.4409e-10  1.042e-11   2.0414e-11 -2.6359e-11; ...
                -3.0881e-12 -2.3602e-12  1.042e-11   8.0458e-13  9.9356e-14 -1.2273e-13; ...
                -4.2085e-13  2.1508e-13  2.0414e-11  9.9356e-14  1.9856e-11 -1.2714e-11; ...
                 1.2413e-11 -4.9476e-12 -2.6359e-11 -1.2273e-13 -1.2714e-11  3.6357e-11]};
            
        seven = {[ 3.748e-11   2.066e-11  -9.1357e-11 -3.0881e-12 -4.2085e-13  1.2413e-11 -3.5535e-08; ...
                   2.066e-11   2.4238e-11 -4.8871e-11 -2.3602e-12  2.1508e-13 -4.9476e-12 -2.2177e-08; ...
                  -9.1357e-11 -4.8871e-11  3.4409e-10  1.042e-11   2.0414e-11 -2.6359e-11  9.2838e-08; ...
                  -3.0881e-12 -2.3602e-12  1.042e-11   8.0458e-13  9.9356e-14 -1.2273e-13  5.4482e-09; ...
                  -4.2085e-13  2.1508e-13  2.0414e-11  9.9356e-14  1.9856e-11 -1.2714e-11 -1.6628e-09; ...
                   1.2413e-11 -4.9476e-12 -2.6359e-11 -1.2273e-13 -1.2714e-11  3.6357e-11 -1.7573e-08; ...
                  -3.5535e-08 -2.2177e-08  9.2838e-08  5.4482e-09 -1.6628e-09 -1.7573e-08  0.00036231]};        
              
        eight = {[ 3.748e-11   2.066e-11  -9.1357e-11 -3.0881e-12 -4.2085e-13  1.2413e-11 -3.5535e-08 0; ...
                   2.066e-11   2.4238e-11 -4.8871e-11 -2.3602e-12  2.1508e-13 -4.9476e-12 -2.2177e-08 0; ...
                  -9.1357e-11 -4.8871e-11  3.4409e-10  1.042e-11   2.0414e-11 -2.6359e-11  9.2838e-08 0; ...
                  -3.0881e-12 -2.3602e-12  1.042e-11   8.0458e-13  9.9356e-14 -1.2273e-13  5.4482e-09 0; ...
                  -4.2085e-13  2.1508e-13  2.0414e-11  9.9356e-14  1.9856e-11 -1.2714e-11 -1.6628e-09 0; ...
                   1.2413e-11 -4.9476e-12 -2.6359e-11 -1.2273e-13 -1.2714e-11  3.6357e-11 -1.7573e-08 0; ...
                  -3.5535e-08 -2.2177e-08  9.2838e-08  5.4482e-09 -1.6628e-09 -1.7573e-08  0.00036231 0; ...
                   0           0           0           0           0           0           0          0]};
            
        nine = {[ 3.748e-11   2.066e-11  -9.1357e-11 -3.0881e-12 -4.2085e-13  1.2413e-11 -3.5535e-08 0  1.3041e-07; ...
                  2.066e-11   2.4238e-11 -4.8871e-11 -2.3602e-12  2.1508e-13 -4.9476e-12 -2.2177e-08 0  8.7873e-08; ...
                 -9.1357e-11 -4.8871e-11  3.4409e-10  1.042e-11   2.0414e-11 -2.6359e-11  9.2838e-08 0 -3.1272e-07; ...
                 -3.0881e-12 -2.3602e-12  1.042e-11   8.0458e-13  9.9356e-14 -1.2273e-13  5.4482e-09 0 -1.2461e-08; ...
                 -4.2085e-13  2.1508e-13  2.0414e-11  9.9356e-14  1.9856e-11 -1.2714e-11 -1.6628e-09 0 -1.9713e-08; ...
                  1.2413e-11 -4.9476e-12 -2.6359e-11 -1.2273e-13 -1.2714e-11  3.6357e-11 -1.7573e-08 0  6.0969e-08; ...
                 -3.5535e-08 -2.2177e-08  9.2838e-08  5.4482e-09 -1.6628e-09 -1.7573e-08  0.00036231 0 -0.00014813; ...
                  0           0           0           0           0           0           0          0  0; ...
                  1.3041e-07  8.7873e-08 -3.1272e-07 -1.2461e-08 -1.9713e-08  6.0969e-08 -0.00014813 0  0.0011749]};
        
        commonOut = {[ 1.75574609146441e-10 -1.76391822498568e-10  1.26701582864723e-11 -8.50798921278566e-12  2.80465474421549e-11  1.78603570386593e-12  5.90834321905617e-08 0 -1.90574918760627e-07; ...
                      -1.76391822498568e-10  2.19847186431344e-10 -2.11714974937397e-11  4.88900720193951e-12 -2.07150417963915e-11  2.63481195571044e-12 -8.18741459439353e-08 0  2.89052907957242e-07; ...
                       1.26701582864723e-11 -2.11714974937397e-11  1.03862044222145e-11 -8.87792456013774e-12  3.61060060754719e-12 -8.98899082284049e-15  1.33873845292676e-08 0 -5.14960908209782e-08; ...
                      -8.50798921278566e-12  4.88900720193951e-12 -8.87792456013774e-12  1.81815864084123e-11 -1.11847510974009e-11 -2.74005163036409e-12 -1.47578050245331e-08 0  3.62925053394288e-08; ...
                       2.80465474421549e-11 -2.07150417963915e-11  3.61060060754719e-12 -1.11847510974009e-11  3.46915046391962e-11  1.057403771701e-11    1.07019379593921e-08 0 -5.41768749409939e-08; ...
                       1.78603570386593e-12  2.63481195571044e-12 -8.98899082284049e-15 -2.74005163036409e-12  1.05740377170100e-11  4.14448895239156e-12 -2.98895146350504e-09 0 -2.96953343396525e-09; ...
                       5.90834321905617e-08 -8.18741459439353e-08  1.33873845292676e-08 -1.47578050245331e-08  1.07019379593921e-08 -2.98895146350504e-09  0.00036231           0 -0.00014813; ...
                       0                     0                     0                     0                     0                     0                     0                    0  0; ...
                      -1.90574918760627e-07  2.89052907957242e-07 -5.14960908209782e-08  3.62925053394288e-08 -5.41768749409939e-08 -2.96953343396525e-09 -0.00014813           0  0.0011749]};
    end

    methods (Test)
        % Tests that 3x3 covariance matrices result in 3x3 outputs
        function test3x3Covariance (testCase, three, defaultR, defaultV, threeOut, relTol)
            RIC = ECI2RIC(three, defaultR, defaultV);
            actSize = size(RIC);
            expSize = [3 3];
            testCase.verifyEqual(actSize, expSize);
            for i = 1:expSize(1)
                for j = 1:expSize(2)
                    testCase.verifyEqual(RIC(i,j),threeOut(i,j),'RelTol',relTol);
                end
            end
        end
        
        % Tests that 6x6 covariance matrices result in 6x6 outputs
        function test6x6Covariance (testCase, six, defaultR, defaultV, commonOut, relTol)
            RIC = ECI2RIC(six, defaultR, defaultV);
            actSize = size(RIC);
            expSize = [6 6];
            testCase.verifyEqual(actSize, expSize);
            for i = 1:expSize(1)
                for j = 1:expSize(2)
                    testCase.verifyEqual(RIC(i,j),commonOut(i,j),'RelTol',relTol);
                end
            end
        end
        
        % Tests that 7x7 covariance matrices result in 7x7 outputs
        function test7x7Covariance (testCase, seven, defaultR, defaultV, commonOut, relTol)
            RIC = ECI2RIC(seven, defaultR, defaultV);
            actSize = size(RIC);
            expSize = [7 7];
            testCase.verifyEqual(actSize, expSize);
            for i = 1:expSize(1)
                for j = 1:expSize(2)
                    testCase.verifyEqual(RIC(i,j),commonOut(i,j),'RelTol',relTol);
                end
            end
        end
        
        % Tests that 8x8 covariance matrices result in 8x8 outputs
        function test8x8Covariance (testCase, eight, defaultR, defaultV, commonOut, relTol)
            RIC = ECI2RIC(eight, defaultR, defaultV);
            actSize = size(RIC);
            expSize = [8 8];
            testCase.verifyEqual(actSize, expSize);
            for i = 1:expSize(1)
                for j = 1:expSize(2)
                    testCase.verifyEqual(RIC(i,j),commonOut(i,j),'RelTol',relTol);
                end
            end
        end
        
        % Tests that 9x9 covariance matrices result in 9x9 outputs
        function test9x9Covariance (testCase, nine, defaultR, defaultV, commonOut, relTol)
            RIC = ECI2RIC(nine, defaultR, defaultV);
            actSize = size(RIC);
            expSize = [9 9];
            testCase.verifyEqual(actSize, expSize);
            for i = 1:expSize(1)
                for j = 1:expSize(2)
                    testCase.verifyEqual(RIC(i,j),commonOut(i,j),'RelTol',relTol);
                end
            end
        end
        
        % Tests that symmetric outputs are found when running in symmetric
        % mode
        function testSymmetry (testCase, nine, defaultR, defaultV)
            RIC = ECI2RIC(nine, defaultR, defaultV);
            nonSymmetryFound = false;
            for i = 1:9
                for j = 1:9
                    if RIC(i,j) ~= RIC(j,i)
                        nonSymmetryFound = true;
                    end
                end
            end
            testCase.verifyEqual(nonSymmetryFound,false);
        end
        
        % Tests that non-symmetric outputs are found when not running in
        % symmetric mode
        function testNonSymmetry (testCase, nine, defaultR, defaultV)
            RIC = ECI2RIC(nine, defaultR, defaultV, false);
            nonSymmetryFound = false;
            for i = 1:9
                for j = 1:9
                    if RIC(i,j) ~= RIC(j,i)
                        nonSymmetryFound = true;
                    end
                end
            end
            testCase.verifyEqual(nonSymmetryFound,true);
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
% E. White       | 08-07-2023 | Initial development
% L. Baars       | 04-23-2025 | Added verification of actual values
%                               produced by the transformations.
% L. Baars       | 10-06-2025 | Updated path fixtures for more consistent
%                               passing runs.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

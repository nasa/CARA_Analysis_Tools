classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        RIC2ECI_UnitTest < matlab.unittest.TestCase
% RIC2ECI - Unit test for RIC2ECI
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Aug 2023;  Latest update: Aug 2023
%
% ----------------- BEGIN CODE -----------------
    
    properties (TestParameter)
        
        defaultR = {[378.39559 4305.721887 5752.767554]};
        defaultV = {[2.360800244 5.580331936 -4.322349039]};

        three = {[ 44.5757544811362   81.6751751052616  -67.8687662707124;
                   81.6751751052616  158.453402956163  -128.616921644857;
                  -67.8687662707124 -128.616921644858   105.490542562701]};
              
        six = {[+3.612023799984759e+03 -9.400236231297700e+03 -6.298653147777720e+02 -1.202684255279736e+00 -1.195413704199584e+00 +1.029005367422424e+01; ...
                -9.400236231297700e+03 +2.450614906178152e+04 +1.607596675826002e+03 +3.121860335933878e+00 +3.148734627832713e+00 -2.679241909910183e+01; ...
                -6.298653147777720e+02 +1.607596675826002e+03 +1.406314287206610e+02 +2.212445034390649e-01 +1.748905714166198e-01 -1.790030136624526e+00; ...
                -1.202684255279736e+00 +3.121860335933878e+00 +2.212445034390649e-01 +4.100943787418265e-04 +3.881475523378179e-04 -3.426808331436135e-03; ...
                -1.195413704199584e+00 +3.148734627832713e+00 +1.748905714166198e-01 +3.881475523378179e-04 +4.338196284756884e-04 -3.412770719926683e-03; ...
                +1.029005367422424e+01 -2.679241909910183e+01 -1.790030136624526e+00 -3.426808331436135e-03 -3.412770719926683e-03 +2.932410138619841e-02]};
            
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
    end

    methods (Test)
        % Tests that 3x3 covariance matrices result in 3x3 outputs
        function test3x3Covariance (testCase, three, defaultR, defaultV)
            act = size(ECI2RIC(three, defaultR, defaultV));
            exp = [3 3];
            testCase.verifyEqual(act, exp);
        end
        
        % Tests that 6x6 covariance matrices result in 6x6 outputs
        function test6x6Covariance (testCase, six, defaultR, defaultV)
            act = size(ECI2RIC(six, defaultR, defaultV));
            exp = [6 6];
            testCase.verifyEqual(act, exp);
        end
        
        % Tests that 7x7 covariance matrices result in 7x7 outputs
        function test7x7Covariance (testCase, seven, defaultR, defaultV)
            act = size(ECI2RIC(seven, defaultR, defaultV));
            exp = [7 7];
            testCase.verifyEqual(act, exp);
        end
        
        % Tests that 8x8 covariance matrices result in 8x8 outputs
        function test8x8Covariance (testCase, eight, defaultR, defaultV)
            act = size(ECI2RIC(eight, defaultR, defaultV));
            exp = [8 8];
            testCase.verifyEqual(act, exp);
        end
        
        % Tests that 9x9 covariance matrices result in 9x9 outputs
        function test9x9Covariance (testCase, nine, defaultR, defaultV)
            act = size(ECI2RIC(nine, defaultR, defaultV));
            exp = [9 9];
            testCase.verifyEqual(act, exp);
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

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

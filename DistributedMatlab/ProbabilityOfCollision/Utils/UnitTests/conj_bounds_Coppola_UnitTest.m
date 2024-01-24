classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        conj_bounds_Coppola_UnitTest < matlab.unittest.TestCase
% conj_bounds_Coppola_UnitTest - Unit test for conj_bounds_Coppola
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
    
    properties (Constant)
        r1 = [-2622716.15829145    1698990.97254698    7053401.88278209];
        v1 = [-4128.50471893557   -5878.79564655151    -118.60127274469];
        C1 = [ 30945.8960599633     44172.992380352      1094.127029015    -16.9384961044964     11.1386598710759     45.9406263498005
                44172.992380352    63081.0350093601    1600.54857635509    -24.1640670703954     15.9361771338151     65.5984249372085
                 1094.127029015    1600.54857635509    110.587616642047   -0.568332213430976    0.456096601496478     1.65805084742508
              -16.9384961044964   -24.1640670703954  -0.568332213430976  0.00928784766563245 -0.00607332603626678  -0.0251311636623764
               11.1386598710759    15.9361771338151   0.456096601496478 -0.00607332603626678  0.00406599298830971   0.0165654694968024
               45.9406263498005    65.5984249372085    1.65805084742508  -0.0251311636623764   0.0165654694968024   0.0682232639995689];
        r2 = [-2631783.12714667    1685565.41148418    7053328.56776558];
        v2 = [-4112.58245687929   -5889.50355090998   -126.881531658755];
        C2 = [ 37176528.2096518    51779295.1953059    1237287.58668701   -20025.8397475978      12497.8621510473      54618.847730921
                51779295.195306    72120313.4426337    1723250.32582261   -27893.1756367764      17408.3723881467     76074.0995332777
               1237287.58668701    1723250.32582261    41984.4552410924   -666.613016373449      416.862049452506      1817.7275694543
              -20025.8397475978   -27893.1756367764   -666.613016373449    10.7880803834256     -6.73318789560721     -29.422096492432
               12497.8621510473    17408.3723881467    416.862049452506   -6.73318789560721      4.20345007734003     18.3622071009157
                54618.847730921    76074.0995332777     1817.7275694543    -29.422096492432      18.3622071009157     80.2453907683188];
    end

    properties (TestParameter)
        defaultR = {(conj_bounds_Coppola_UnitTest.r2 - conj_bounds_Coppola_UnitTest.r1)'};
        defaultV = {(conj_bounds_Coppola_UnitTest.v2 - conj_bounds_Coppola_UnitTest.v1)'};
        defaultC = {conj_bounds_Coppola_UnitTest.C1 + conj_bounds_Coppola_UnitTest.C2};
        defaultHBR = {14.8};
        defaultGamma = {1e-16};
    end

    methods (Test)
        % Tests that output is the same without an argument and with its
        % default value
        function testDefaultParameter (testCase, defaultGamma, defaultHBR, defaultR, defaultV, defaultC)
            [a_t0, a_t1, a_t0g, a_t1g] = conj_bounds_Coppola(defaultGamma, defaultHBR, defaultR, defaultV, defaultC);
            [e_t0, e_t1, e_t0g, e_t1g] = conj_bounds_Coppola(defaultGamma, defaultHBR, defaultR, defaultV, defaultC, false);
            
            act = [a_t0, a_t1, a_t0g, a_t1g];
            exp = [e_t0, e_t1, e_t0g, e_t1g];
            
            testCase.verifyEqual(act, exp);
        end
        
        % Tests that passing vector inputs returns outputs of the correct
        % dimensions
        function testVectorizedOutputDimensions (testCase, defaultGamma, defaultHBR, defaultR, defaultV, defaultC)
            gamma = repmat(defaultGamma, 3, 2);
            [t0, t1, t0g, t1g] = conj_bounds_Coppola(gamma, defaultHBR, defaultR, defaultV, defaultC);
            
            testCase.verifyEqual(size(t0), size(gamma));
            testCase.verifyEqual(size(t1), size(gamma));
            testCase.verifyEqual(size(t0g), [1 1]);
            testCase.verifyEqual(size(t1g), [1 1]);
        end
        
        % Tests that both 3x3 and 6x6 covariance matrices can be used
        % without error (note that since these Pc calculations assume no
        % velocity covariance, the cases will be equal)
        function testCovarianceSizes (testCase, defaultGamma, defaultHBR, defaultR, defaultV, defaultC)
            C3 = defaultC(1:3, 1:3);
            [t03, t13, t0g3, t1g3] = conj_bounds_Coppola(defaultGamma, defaultHBR, defaultR, defaultV, C3);
            [t06, t16, t0g6, t1g6] = conj_bounds_Coppola(defaultGamma, defaultHBR, defaultR, defaultV, defaultC);
            
            act = [t03, t13, t0g3, t1g3];
            exp = [t06, t16, t0g6, t1g6];
            
            testCase.verifyEqual(act, exp);
        end
        
        % Tests bounds for zero relative velocity
        function testZeroRelativeVelocity (testCase, defaultGamma, defaultHBR, defaultR, defaultC)
            [t0, t1, t0g, t1g] = conj_bounds_Coppola(defaultGamma, defaultHBR, defaultR, [0, 0, 0], defaultC);
            act = [t0, t1, t0g, t1g];
            exp = [-Inf, Inf, -Inf, Inf];
            
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

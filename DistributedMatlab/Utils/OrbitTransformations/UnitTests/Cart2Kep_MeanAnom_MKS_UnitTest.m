classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Cart2Kep_MeanAnom_MKS_UnitTest < matlab.unittest.TestCase
% Cart2Kep_MeanAnom_MKS_UnitTest - Unit test for Cart2Kep_MeanAnom_MKS
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Jun 2023;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------   
    
    properties (Constant)
        mu = 3.986004418e14;
    end

    properties (TestParameter)
        % Verification case from David Vallado's Fundamentals of Astrodynamics and Applications
        r = {[6524834, 6862875, 6448296]};
        v = {[4901.327, 5533.756, -1976.341]};
        
        res = {[36127337, 0.832853, 1.5336, 3.9774, 0.9317, 0.13273]};
        
        zeroVersR = { ...
            [0, 6600000, 6700000], ...
            [6500000, 0, 6700000], ...
            [6500000, 6600000, 0], ...
            [6500000, 6600000, 6700000], ...
            [6500000, 6600000, 6700000], ...
            [6500000, 6600000, 6700000] ...
        };
    
        zeroVersV = { ...
            [1000, 2000, 3000], ...
            [1000, 2000, 3000], ...
            [1000, 2000, 3000], ...
            [0, 2000, 3000], ...
            [1000, 0, 3000], ...
            [1000, 2000, 0] ...
        }
    
        zeroResults = {...
            [5632700.001, 0.970753, 2.189742, 2.379091, 4.310613, 1.629306], ...
            [5582623.190, 0.877344, 0.959048, 5.474890, 4.396141, 1.879656], ...
            [5531547.202, 0.794455, 1.344445, 0.793032, 3.334298, 2.165707], ...
            [7026069.266, 0.888400, 1.006163, 0.317128, 4.109434, 1.725277], ...
            [6673187.117, 0.880137, 1.843742, 0.996910, 3.958966, 1.999724], ...
            [6157736.446, 0.944387, 1.167073, 4.248741, 5.666857, 2.309350] ...
        };
    end
    
    methods (Test)
        % Tests against known case from Vallado
        function testKnownCase (testCase, r, v, res)
            [a, e, i, Omega, w, M] = Cart2Kep_MeanAnom_MKS(r, v);
            testCase.verifyEqual([a, e, i, Omega, w, M], res, 'RelTol', 1E-4, 'AbsTol', 1E-4);
        end
    end
    
    methods (Test, ParameterCombination = 'sequential')
        % Tests zeroes in Cartesian elements input
        function testCartesianElementZeroes (testCase, zeroVersR, zeroVersV, zeroResults)
            [a, e, i, Omega, w, M] = Cart2Kep_MeanAnom_MKS(zeroVersR, zeroVersV);
            testCase.verifyEqual([a, e, i, Omega, w, M], zeroResults, 'RelTol', 1E-5, 'AbsTol', 1E-5);
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
% E. White       | 06-12-2023 | Initial development
% E. White       | 06-13-2023 | Fixed error in documentation
% E. White       | 06-30-2023 | Added AbsTol check for zero values

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

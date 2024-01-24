classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Kep2Cart_MeanAnom_MKS_UnitTest < matlab.unittest.TestCase
% Kep2Cart_MeanAnom_MKS_UnitTest - Unit test for Kep2Cart_MeanAnom_MKS
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
        ver = {[36126642.83, 0.83285, 1.53362081372742, 3.97743083236988, 0.931656754714573, 0.132731244829756]};
        resR = {[6525368.120986; 6861531.834896; 6449118.614160]};
        resV = {[4902.2786446; 5533.1395663; -1975.7100988]};
  
        specialVers = { ...
            [6500000, 0, 0, 0, 0, 0], ... % Equatorial circular, r_j > 0
            [6500000, 0, 0, 0, 0, 3 * pi / 2], ... % Equatorial circular, r_j < 0
            [6500000, 0, pi / 4, 0, 0, 0], ... % Inclined circular, n_j > 0, r_k > 0
            [4082084.867, 0.5923235828, 0, 0, pi, pi] ... % Equatorial elliptical, e_j > 0, r dot v > 0
        };
    
        specialResR = { ...
            [6500000, 0, 0], ... % Equatorial circular, r_j > 0
            [0, -6500000, 0], ... % Equatorial circular, r_j < 0
            [6500000, 0, 0], ... % Inclined circular, n_j > 0, r_k > 0
            [6500000, 0, 0] ... % Equatorial elliptical, e_j > 0, r dot v > 0
        };
    
            specialResV = { ...
            [0, sqrt(Kep2Cart_MeanAnom_MKS_UnitTest.mu / 6500000), 0], ... % Equatorial circular, r_j > 0
            [sqrt(Kep2Cart_MeanAnom_MKS_UnitTest.mu / 6500000), 0, 0], ... % Equatorial circular, r_j < 0
            [0, sqrt(Kep2Cart_MeanAnom_MKS_UnitTest.mu / 13000000), sqrt(Kep2Cart_MeanAnom_MKS_UnitTest.mu / 13000000)], ... % Inclined circular, n_j > 0, r_k > 0
            [0, 5000, 0] ... % Equatorial elliptical, e_j > 0, r dot v > 0
        };
    end
    
    methods (Test)
        % Tests against known case from Vallado
        function testKnownCase (testCase, ver, resR, resV)
            a = ver(1);
            e = ver(2);
            i = ver(3);
            Omega = ver(4);
            w = ver(5);
            M = ver(6);
            
            [r, v] = Kep2Cart_MeanAnom_MKS(a, e, i, Omega, w, M);
            testCase.verifyEqual([r, v], [resR, resV], 'RelTol', 1E-4, 'AbsTol', 1E-4);
        end
    end
    
    methods (Test, ParameterCombination = 'sequential')
        % Tests special case orbits
        function testSpecialCaseOrbits (testCase, specialVers, specialResR, specialResV)
            a = specialVers(1);
            e = specialVers(2);
            i = specialVers(3);
            Omega = specialVers(4);
            w = specialVers(5);
            M = specialVers(6);
            
            [r, v] = Kep2Cart_MeanAnom_MKS(a, e, i, Omega, w, M);
            testCase.verifyEqual([r, v], [specialResR', specialResV'], 'RelTol', 1E-4, 'AbsTol', 1E-4);
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
% E. White       | 06-13-2023 | Initial development

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

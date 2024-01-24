classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Eq2Kep_UnitTest < matlab.unittest.TestCase
% Eq2Kep_UnitTest - Unit test for Eq2Kep
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
        mu = 3.986004418e5;
        maxR = 43000;
        minR = 6400;
        defaultA = 6500;
    end
    
    properties (TestParameter)
        ver = {[-0.186939977644424, 0.0251358447064123, 5.19960068977612, 0.000618165872034742, -0.206130724502682, -2.36641891694603]};
        verRes = {[10141.6645636734, 0.188622283759906, 2.344662081055977,3.228479950052446,6.062640205399887,2.191665841503368]};
        
        parabolic = {[1, 0, 0, sqrt(Eq2Kep_UnitTest.mu / Eq2Kep_UnitTest.defaultA ^ 3), 0, 0]};
        rectilinear = {[1, 0, 0, 0, 0, 0]};
        
        equatorialRetrograde = {[0, 0, 0, sqrt(Eq2Kep_UnitTest.mu / Eq2Kep_UnitTest.defaultA ^ 3), 0, 0, -1]};
        equatorialRetrogradeVer = {[Eq2Kep_UnitTest.defaultA, 0, pi, 0, 0, 0]};
    end

    methods (Test)
        % Tests parabolic and rectilinear orbits
        function testParabolicAndRectilinearOrbits (testCase, parabolic, rectilinear)
            testCase.verifyError(@() Eq2Kep(parabolic), 'Eq2Kep:ParabolicOrRectilinearOrbit');
            testCase.verifyError(@() Eq2Kep(rectilinear), 'Eq2Kep:ParabolicOrRectilinearOrbit');
        end

        % TODO: Add hyperbolic orbit test
        
        % Tests equatorial retrograde orbits
        function testEquatorialRetrogradeOrbit (testCase, equatorialRetrograde, equatorialRetrogradeVer)
            testCase.verifyEqual(Eq2Kep(equatorialRetrograde), equatorialRetrogradeVer, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests known verification case
        function testKnownCase (testCase, ver, verRes)
            testCase.verifyEqual(Eq2Kep(ver), verRes, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests Equinoctial > Keplerian > Equinoctial > Keplerian >
        % Equinoctial
        function testTwoWay (testCase)
            for i = 1:1000
                r_1 = (Eq2Kep_UnitTest.maxR - Eq2Kep_UnitTest.minR) * rand + Eq2Kep_UnitTest.minR;
                r_2 = (Eq2Kep_UnitTest.maxR - Eq2Kep_UnitTest.minR) * rand + Eq2Kep_UnitTest.minR;

                a = (r_1 + r_2) / 2;
                e = rand;

                r_p = a * (1 - e);

                while r_p < Eq2Kep_UnitTest.minR
                    e = rand;
                    r_p = a * (1 - e);
                end

                vec = [a, e, pi * rand, 2 * pi * rand, 2 * pi * rand, 2 * pi * rand];
                eq = Kep2Eq(vec);

                testCase.verifyEqual(Kep2Eq(Eq2Kep(Kep2Eq(Eq2Kep(eq)))), eq, 'RelTol', 1E-5, 'AbsTol', 1E-5);
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
% E. White       | 06-30-2023 | Initial development

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

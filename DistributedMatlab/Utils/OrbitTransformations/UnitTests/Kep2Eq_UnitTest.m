classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Kep2Eq_UnitTest < matlab.unittest.TestCase
% Kep2Eq_UnitTest - Unit test for Kep2Eq
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
        ver = {[10141.6645636734, 0.188622283759906, 2.344662081055977,3.228479950052446,6.062640205399887,2.191665841503368]};
        verRes = {[-0.186939977644424, 0.0251358447064123, 5.19960068977612, 0.000618165872034742, -0.206130724502682, -2.36641891694603, 1]};
        
        parabolic = {[inf, 1, 0, 0, 0, 0]};
        rectilinear = {[Kep2Eq_UnitTest.defaultA, 1, 0, 0, 0, 0]};
        hyperbolic = {[-Kep2Eq_UnitTest.defaultA, 1.5, 0, 0, 0, 0]};
        
        equatorialRetrograde = {[Kep2Eq_UnitTest.defaultA, 0, pi, 0, 0, 0]};
        equatorialRetrogradeVer = {[0, 0, 0, sqrt(Kep2Eq_UnitTest.mu / Kep2Eq_UnitTest.defaultA ^ 3), 0, 0, -1]};
    end

    methods (Test)
        % Tests parabolic and rectilinear orbits
        function testParabolicAndRectilinearOrbits (testCase, parabolic, rectilinear)
            testCase.verifyError(@() Kep2Eq(parabolic), 'Kep2Eq:ParabolicOrRectilinearOrbit');
            testCase.verifyError(@() Kep2Eq(rectilinear), 'Kep2Eq:ParabolicOrRectilinearOrbit');
        end
        
        % Tests hyperbolic orbits
        function testHyperbolicOrbit (testCase, hyperbolic)
            testCase.verifyError(@() Kep2Eq(hyperbolic), 'Kep2Eq:HyperbolicOrbit');
        end
        
        % Tests equatorial retrograde orbits
        function testEquatorialRetrogradeOrbit (testCase, equatorialRetrograde, equatorialRetrogradeVer)
            testCase.verifyEqual(Kep2Eq(equatorialRetrograde), equatorialRetrogradeVer, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests known verification case
        function testKnownCase (testCase, ver, verRes)
            testCase.verifyEqual(Kep2Eq(ver), verRes, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests Keplerian > Equinoctial > Keplerian > Equinoctial >
        % Keplerian
        function testTwoWay (testCase)
            for i = 1:1000
                r_1 = (Kep2Eq_UnitTest.maxR - Kep2Eq_UnitTest.minR) * rand + Kep2Eq_UnitTest.minR;
                r_2 = (Kep2Eq_UnitTest.maxR - Kep2Eq_UnitTest.minR) * rand + Kep2Eq_UnitTest.minR;

                a = (r_1 + r_2) / 2;
                e = rand;

                r_p = a * (1 - e);

                while r_p < Kep2Eq_UnitTest.minR
                    e = rand;
                    r_p = a * (1 - e);
                end

                vec = [a, e, pi * rand, 2 * pi * rand, 2 * pi * rand, 2 * pi * rand];

                testCase.verifyEqual(Eq2Kep(Kep2Eq(Eq2Kep(Kep2Eq(vec)))), vec, 'RelTol', 1E-5, 'AbsTol', 1E-5);
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

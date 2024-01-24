classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        convert_cartesian_to_equinoctial_UnitTest < matlab.unittest.TestCase
% convert_cartesian_to_equinoctial_UnitTest - Unit test for convert_cartesian_to_equinoctial
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
        defaultR = 6500;
    end
    
    properties (TestParameter)
        ver = {[6513.711946044443, 6882.748803674651, 6438.350633995484, 4.902960531245754, -1.518827797965108, -1.983803852683226]};
        verRes = {[10141.66456, 0.000618165872034742, -0.186939977644424, 0.0251358447064123, -0.206130724502682, -2.36641891694603, -1.08358461740347]};
        
        parabolic = {[convert_cartesian_to_equinoctial_UnitTest.defaultR, 0, 0, 0, sqrt(2 * convert_cartesian_to_equinoctial_UnitTest.mu / convert_cartesian_to_equinoctial_UnitTest.defaultR), 0]};
        hyperbolic = {[convert_cartesian_to_equinoctial_UnitTest.defaultR, 0, 0, 0, 2 * sqrt(2 * convert_cartesian_to_equinoctial_UnitTest.mu / convert_cartesian_to_equinoctial_UnitTest.defaultR), 0]};
        
        equatorial = {[convert_cartesian_to_equinoctial_UnitTest.defaultR, 0, 0, 0, sqrt(convert_cartesian_to_equinoctial_UnitTest.mu / convert_cartesian_to_equinoctial_UnitTest.defaultR), 0]};
        equatorialRetrograde = {[convert_cartesian_to_equinoctial_UnitTest.defaultR, 0, 0, 0, -sqrt(convert_cartesian_to_equinoctial_UnitTest.mu / convert_cartesian_to_equinoctial_UnitTest.defaultR), 0]};
    end

    methods (Test)
        % Tests default parameters
        function testDefaultParameters (testCase, ver)
            [actA,actN,actAf,actAg,actChi,actPsi,actLM,actF] = convert_cartesian_to_equinoctial(ver(1:3), ver(4:6));
            [expA,expN,expAf,expAg,expChi,expPsi,expLM,expF] = convert_cartesian_to_equinoctial(ver(1:3), ver(4:6), 1, convert_cartesian_to_equinoctial_UnitTest.mu, true);
            act = [actA,actN,actAf,actAg,actChi,actPsi,actLM,actF];
            exp = [expA,expN,expAf,expAg,expChi,expPsi,expLM,expF];
            testCase.verifyEqual(act, exp);
        end

        % Tests invalid retrograde factor
        function testInvalidRetrogradeFactor (testCase, ver)
            warning('on');
            testCase.verifyWarning(@() convert_cartesian_to_equinoctial(ver(1:3), ver(4:6), 2), 'convert_cartesian_to_equinoctial:InvalidRetrogradeFactor');
        end

        % Tests unbound orbits
        function testUnboundOrbits (testCase, parabolic, hyperbolic)
            warning('on');
            testCase.verifyWarning(@() convert_cartesian_to_equinoctial(parabolic(1:3), parabolic(4:6)), 'convert_cartesian_to_equinoctial:UnboundOrbit');
            testCase.verifyWarning(@() convert_cartesian_to_equinoctial(hyperbolic(1:3), hyperbolic(4:6)), 'convert_cartesian_to_equinoctial:UnboundOrbit');
        end

        % Tests singular case (i = 180 deg with fr = -1 or i = 0 deg with fr = -1)
        function testSingularity (testCase, equatorialRetrograde, equatorial)
            warning('on');
            testCase.verifyWarning(@() convert_cartesian_to_equinoctial(equatorialRetrograde(1:3), equatorialRetrograde(4:6)), 'convert_cartesian_to_equinoctial:EquatorialRetrogradeOrbit');
            testCase.verifyWarning(@() convert_cartesian_to_equinoctial(equatorial(1:3), equatorial(4:6), -1), 'convert_cartesian_to_equinoctial:EquatorialRetrogradeOrbit');
        end
        
        % Tests known verification case
        
        % Ignores final return value (F) as it is an intermediate
        % calculation and other values would be incorrect if it were
        % incorrect
        function testKnownCase (testCase, ver, verRes)
            [a,n,af,ag,chi,psi,lM] = convert_cartesian_to_equinoctial(ver(1:3), ver(4:6));
            testCase.verifyEqual([a,n,af,ag,chi,psi,lM], verRes, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests Cartesian > Equinoctial > Cartesian > Equinoctial >
        % Cartesian
        function testTwoWay (testCase)
            for i = 1:1000
                r_1 = (convert_cartesian_to_equinoctial_UnitTest.maxR - convert_cartesian_to_equinoctial_UnitTest.minR) * rand + convert_cartesian_to_equinoctial_UnitTest.minR;
                r_2 = (convert_cartesian_to_equinoctial_UnitTest.maxR - convert_cartesian_to_equinoctial_UnitTest.minR) * rand + convert_cartesian_to_equinoctial_UnitTest.minR;

                a = (r_1 + r_2) / 2;
                e = rand;

                r_p = a * (1 - e);

                while r_p < convert_cartesian_to_equinoctial_UnitTest.minR
                    e = rand;
                    r_p = a * (1 - e);
                end

                vec = [a, e, pi * rand, 2 * pi * rand, 2 * pi * rand, 2 * pi * rand];
                cart = Kep2Cart(vec);
                [~,n,af,ag,chi,psi,lM] = convert_cartesian_to_equinoctial(cart(1:3), cart(4:6));
                [r, v] = convert_equinoctial_to_cartesian(n,af,ag,chi,psi,lM, 0);
                [~,n,af,ag,chi,psi,lM] = convert_cartesian_to_equinoctial(r, v);
                [r, v] = convert_equinoctial_to_cartesian(n,af,ag,chi,psi,lM, 0);

                testCase.verifyEqual([r; v]', cart, 'RelTol', 1E-5, 'AbsTol', 1E-5);
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

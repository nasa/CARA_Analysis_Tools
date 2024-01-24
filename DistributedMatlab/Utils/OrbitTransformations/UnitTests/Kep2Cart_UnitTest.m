classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Kep2Cart_UnitTest < matlab.unittest.TestCase
% Kep2Cart_UnitTest - Unit test for Kep2Cart
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
        defaultR = 6500;
    end

    properties (TestParameter)
        ver = {[36126.64283, 0.83285, 1.53362081372742, 3.97743083236988, 0.931656754714573, 0.132731244829756]};
        res = {[6525.368120986, 6861.531834896, 6449.118614160, 4.9022786446, 5.5331395663, -1.9757100988]};
        
        parabolic = {[inf, 1, 0, 0, 0, 0]};
        hyperbolic = {[Kep2Cart_UnitTest.defaultR, 1.5, 0, 0, 0, 0]};
        rectilinear = {[Kep2Cart_UnitTest.defaultR, 1, 0, 0, 0, 0]};
        
        circular = {[Kep2Cart_UnitTest.defaultR, 0, pi / 4, 0, 0, 0]};
        equatorial = {[4082.084867, 0.5923235828, 0, 0, pi, pi]};
        
        specialVers = { ...
            [Kep2Cart_UnitTest.defaultR, 0, 0, 0, 0, 0], ... % Equatorial circular, r_j >= 0
            [Kep2Cart_UnitTest.defaultR, 0, 0, 0, 0, 3 * pi / 2], ... % Equatorial circular, r_j < 0
            [Kep2Cart_UnitTest.defaultR, 0, pi / 4, 0, 0, 0], ... % Inclined circular, n_j >= 0, r_k >= 0
            [Kep2Cart_UnitTest.defaultR, 0, pi / 4, 3 * pi / 2, 0, 0], ... % Inclined circular, n_j < 0, r_k >= 0
            [Kep2Cart_UnitTest.defaultR, 0, pi / 4, 0, 0, 4.71238898038469], ... % Inclined circular, n_j >= 0, r_k < 0
            [Kep2Cart_UnitTest.defaultR, 0, pi / 4, 3 * pi / 2, 0, 4.71238898038469], ... % Inclined circular, n_j < 0, r_k < 0
            [Kep2Cart_UnitTest.defaultR, sqrt(2) / 2, 0, 0, pi / 4, 0.863689545608349] ... % Equatorial elliptical, e_j >= 0, r . v >= 0
            [Kep2Cart_UnitTest.defaultR, sqrt(2) / 2, 0, 0, 7 * pi / 4, 0.863689545608349] ... % Equatorial elliptical, e_j < 0, r . v >= 0
            [Kep2Cart_UnitTest.defaultR, sqrt(2) / 2, 0, 0, 3 * pi / 4, 5.41949576157124] ... % Equatorial elliptical, e_j >= 0, r . v < 0
            [Kep2Cart_UnitTest.defaultR, sqrt(2) / 2, 0, 0, 5 * pi / 4, 5.41949576157124] ... % Equatorial elliptical, e_j < 0, r . v < 0
        };
    
        specialResults = { ...
            [Kep2Cart_UnitTest.defaultR, 0, 0, 0, sqrt(Kep2Cart_UnitTest.mu / Kep2Cart_UnitTest.defaultR), 0], ... % Equatorial circular, r_j >= 0
            [0, -Kep2Cart_UnitTest.defaultR, 0, sqrt(Kep2Cart_UnitTest.mu / Kep2Cart_UnitTest.defaultR), 0, 0], ... % Equatorial circular, r_j < 0
            [Kep2Cart_UnitTest.defaultR, 0, 0, 0, sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR))], ... % Inclined circular, n_j >= 0, r_k >= 0
            [0, -Kep2Cart_UnitTest.defaultR, 0, sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), 0, sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR))], ... % Inclined circular, n_j < 0, r_k >= 0
            [0, - Kep2Cart_UnitTest.defaultR / sqrt(2), - Kep2Cart_UnitTest.defaultR / sqrt(2), sqrt(Kep2Cart_UnitTest.mu / Kep2Cart_UnitTest.defaultR), 0, 0], ... % Inclined circular, n_j >= 0, r_k < 0
            [-Kep2Cart_UnitTest.defaultR / sqrt(2), 0, - Kep2Cart_UnitTest.defaultR / sqrt(2), 0, -sqrt(Kep2Cart_UnitTest.mu / Kep2Cart_UnitTest.defaultR), 0], ... % Inclined circular, n_j < 0, r_k < 0
            [-Kep2Cart_UnitTest.defaultR, 0, 0, -sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), -sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j >= 0, r . v >= 0
            [0, Kep2Cart_UnitTest.defaultR, 0, -sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j < 0, r . v >= 0
            [Kep2Cart_UnitTest.defaultR, 0, 0, -sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j >= 0, r . v < 0
            [0, Kep2Cart_UnitTest.defaultR, 0, -sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), -sqrt(Kep2Cart_UnitTest.mu / (2 * Kep2Cart_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j < 0, r . v < 0
        };
    
        arrVer = {...
            [ ...
                36126.64283, 0.83285, 1.53362081372742, 3.97743083236988, 0.931656754714573, 0.132731244829756; ...
                4082.084867, 0.5923235828, 0, 0, pi, pi ...
            ] ...
        };
    
        arrResults = { ...
            [ ...
                6525.368120986, 6861.531834896, 6449.118614160, 4.9022786446, 5.5331395663, -1.9757100988; ...
                6500, 0, 0, 0, 5, 0 ...
            ] ...
        };
    
    end
    
    methods (Test)
        % Tests against invalid orbital elements
        function testInvalidOrbitalElements (testCase, ver)
            testCase.verifyError(@() Kep2Cart([]), 'Kep2Cart:InvalidOrbitalElements'); % Too few Keplerian elements
            testCase.verifyError(@() Kep2Cart([6500, 0, 0, 0, 0, 0, 0]), 'Kep2Cart:InvalidOrbitalElements'); % Too many Keplerian elements
            testCase.verifyError(@() Kep2Cart(['a', 'e', 'i', 'Omega', 'w', 'M']), 'Kep2Cart:InvalidOrbitalElements'); % Invalid Keplerian elements
            testCase.verifyError(@() Kep2Cart(ver'), 'Kep2Cart:InvalidOrbitalElements'); % Incorrect format
        end
        
        % Tests against parabolic and rectilinear orbits
        function testParabolicAndRectilinearOrbits (testCase, parabolic, rectilinear)
            testCase.verifyError(@() Kep2Cart(parabolic), 'Kep2Cart:ParabolicOrRectilinearOrbit');
            testCase.verifyError(@() Kep2Cart(rectilinear), 'Kep2Cart:ParabolicOrRectilinearOrbit');
        end
        
        % Tests against hyperbolic orbits
        function testHyperbolicOrbit (testCase, hyperbolic)
            testCase.verifyError(@() Kep2Cart(hyperbolic), 'Kep2Cart:HyperbolicOrbit');
        end
        
        % Tests against known verification case from David Vallado's Fundamentals of Astrodynamics and Applications
        function testKnownCase (testCase, ver, res)
            testCase.verifyEqual(Kep2Cart(ver), res, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
       % Tests warnings for special case orbits
       function testSpecialCaseOrbitWarnings (testCase, circular, equatorial)
           warning('on');
           testCase.verifyWarning(@() Kep2Cart(circular), 'Kep2Cart:SpecialCaseOrbit'); % Circular orbit
           testCase.verifyWarning(@() Kep2Cart(equatorial), 'Kep2Cart:SpecialCaseOrbit'); % Equatorial orbit
       end
       
       % Tests array inputs
       function testArrayInputs (testCase, arrVer, arrResults)
           testCase.verifyEqual(Kep2Cart(arrVer), arrResults, 'RelTol', 1E-5, 'AbsTol', 1E-5);
       end
       
        % Tests Cartesian > Keplerian > Cartesian > Keplerian > Cartesian
        function testTwoWay (testCase)
            for i = 1:1000
                r_1 = (Cart2Kep_UnitTest.maxR - Cart2Kep_UnitTest.minR) * rand + Cart2Kep_UnitTest.minR;
                r_2 = (Cart2Kep_UnitTest.maxR - Cart2Kep_UnitTest.minR) * rand + Cart2Kep_UnitTest.minR;

                a = (r_1 + r_2) / 2;
                e = rand;

                r_p = a * (1 - e);

                while r_p < Cart2Kep_UnitTest.minR
                    e = rand;
                    r_p = a * (1 - e);
                end

                vec = [a, e, pi * rand, 2 * pi * rand, 2 * pi * rand, 2 * pi * rand];
                cart = Kep2Cart(vec);
                out = Kep2Cart(Cart2Kep(Kep2Cart(Cart2Kep(cart, 'Mean', 'Rad')), 'Mean', 'Rad'));
                
                testCase.verifyEqual(out, cart, 'RelTol', 1E-5, 'AbsTol', 1E-5);
            end
        end
       
    end
    
    methods (Test, ParameterCombination = 'sequential')
        % Tests special case orbits
        function testSpecialCaseOrbits(testCase, specialVers, specialResults)
            testCase.verifyEqual(Kep2Cart(specialVers), specialResults, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests Cartesian > Keplerian > Cartesian for special cases
        function testSpecialTwoWay (testCase, specialVers)
            testCase.verifyEqual(Cart2Kep(Kep2Cart(specialVers), 'Mean', 'Rad'), specialVers, 'RelTol', 1E-5, 'AbsTol', 1E-5);
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
% E. White       | 06-30-2023 | Added unit tests for remaining
%                               special-case orbits, added test for
%                               hyperbolic orbits, added back and forth
%                               checks for special cases, added randomized
%                               back and forth checks

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

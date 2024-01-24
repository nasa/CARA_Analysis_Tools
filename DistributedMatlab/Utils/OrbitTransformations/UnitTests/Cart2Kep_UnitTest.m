classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Cart2Kep_UnitTest < matlab.unittest.TestCase
% Cart2Kep_UnitTest - Unit test for Cart2Kep
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
        maxR = 43000;
        minR = 6400;
    end

    properties (TestParameter)
        anomTypes = {'True', 'Eccentric', 'Mean'};
        
        % Verification case from David Vallado's Fundamentals of Astrodynamics and Applications
        ver = {[6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341]};
        
        parabolic = {[Cart2Kep_UnitTest.defaultR, 0, 0, 0, sqrt(2 * Cart2Kep_UnitTest.mu / Cart2Kep_UnitTest.defaultR), 0]};
        hyperbolic = {[Cart2Kep_UnitTest.defaultR, 0, 0, 0, 15, 0]};
        rectilinear = {[Cart2Kep_UnitTest.defaultR, 0, 0, 5, 0, 0]};
        
        equatorialCircular = {[Cart2Kep_UnitTest.defaultR, 0, 0, 0, sqrt(Cart2Kep_UnitTest.mu / Cart2Kep_UnitTest.defaultR), 0]};
        inclinedCircular = {[Cart2Kep_UnitTest.defaultR, 0, 0, 0, sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR))]};
        equatorialElliptical = {[Cart2Kep_UnitTest.defaultR, 0, 0, 0, 5, 0]};
        
        degResVec = { ...
            [36127.337, 0.832853, 87.870, 227.89, 53.38, 92.335], ...
            [36127.337, 0.832853, 87.870, 227.89, 53.38, 34.9218704046551], ...
            [36127.337, 0.832853, 87.870, 227.89, 53.38, 7.60471328476942] ...
        };
    
        radResVec = { ...
            [36127.337, 0.832853, 1.5336, 3.9774, 0.9317, 1.61155], ...
            [36127.337, 0.832853, 1.5336, 3.9774, 0.9317, .60950], ...
            [36127.337, 0.832853, 1.5336, 3.9774, 0.9317, 0.13273] ...
        };
        
        zeroVers = { ...
            [0, 6600, 6700, 1, 2, 3], ...
            [6500, 0, 6700, 1, 2, 3], ...
            [6500, 6600, 0, 1, 2, 3], ...
            [6500, 6600, 6700, 0, 2, 3], ...
            [6500, 6600, 6700, 1, 0, 3], ...
            [6500, 6600, 6700, 1, 2, 0] ...
        };
    
        zeroResults = {...
            [5632.700001, 0.970753, 125.462947, 136.311888, 246.979929, 174.024926], ...
            [5582.623190, 0.877344, 54.949377, 313.688112, 251.880348, 169.371138], ...
            [5531.547202, 0.794455, 77.031032, 45.437364, 191.041214, 168.958786], ...
            [7026.069266, 0.888400, 57.648879, 18.170101, 235.453261, 168.472901], ...
            [6673.187117, 0.880137, 105.638663, 57.118760, 226.832062, 170.655005], ...
            [6157.736446, 0.944387, 66.868354, 243.43495, 324.686968, 175.723133] ...
        };
    
        specialVers = { ...
            [Cart2Kep_UnitTest.defaultR, 0, 0, 0, sqrt(Cart2Kep_UnitTest.mu / Cart2Kep_UnitTest.defaultR), 0], ... % Equatorial circular, r_j >= 0
            [0, -Cart2Kep_UnitTest.defaultR, 0, sqrt(Cart2Kep_UnitTest.mu / Cart2Kep_UnitTest.defaultR), 0, 0], ... % Equatorial circular, r_j < 0
            [Cart2Kep_UnitTest.defaultR, 0, 0, 0, sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR))], ... % Inclined circular, n_j >= 0, r_k >= 0
            [0, -Cart2Kep_UnitTest.defaultR, 0, sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), 0, sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR))], ... % Inclined circular, n_j < 0, r_k >= 0
            [0, - Cart2Kep_UnitTest.defaultR / sqrt(2), - Cart2Kep_UnitTest.defaultR / sqrt(2), sqrt(Cart2Kep_UnitTest.mu / Cart2Kep_UnitTest.defaultR), 0, 0], ... % Inclined circular, n_j >= 0, r_k < 0
            [-Cart2Kep_UnitTest.defaultR / sqrt(2), 0, - Cart2Kep_UnitTest.defaultR / sqrt(2), 0, -sqrt(Cart2Kep_UnitTest.mu / Cart2Kep_UnitTest.defaultR), 0], ... % Inclined circular, n_j < 0, r_k < 0
            [-Cart2Kep_UnitTest.defaultR, 0, 0, -sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), -sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j >= 0, r . v >= 0
            [0, Cart2Kep_UnitTest.defaultR, 0, -sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j < 0, r . v >= 0
            [Cart2Kep_UnitTest.defaultR, 0, 0, -sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j >= 0, r . v < 0
            [0, Cart2Kep_UnitTest.defaultR, 0, -sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), -sqrt(Cart2Kep_UnitTest.mu / (2 * Cart2Kep_UnitTest.defaultR)), 0] ... % Equatorial elliptical, e_j < 0, r . v < 0
        };
        
    
        specialResults = { ...
            [Cart2Kep_UnitTest.defaultR, 0, 0, 0, 0, 0], ... % Equatorial circular, r_j >= 0
            [Cart2Kep_UnitTest.defaultR, 0, 0, 0, 0, 270], ... % Equatorial circular, r_j < 0
            [Cart2Kep_UnitTest.defaultR, 0, 45, 0, 0, 0], ... % Inclined circular, n_j >= 0, r_k >= 0
            [Cart2Kep_UnitTest.defaultR, 0, 45, 270, 0, 0], ... % Inclined circular, n_j < 0, r_k >= 0
            [Cart2Kep_UnitTest.defaultR, 0, 45, 0, 0, 270], ... % Inclined circular, n_j >= 0, r_k < 0
            [Cart2Kep_UnitTest.defaultR, 0, 45, 270, 0, 270], ... % Inclined circular, n_j < 0, r_k < 0
            [Cart2Kep_UnitTest.defaultR, sqrt(2) / 2, 0, 0, 45, 135] ... % Equatorial elliptical, e_j >= 0, r . v >= 0
            [Cart2Kep_UnitTest.defaultR, sqrt(2) / 2, 0, 0, 315, 135] ... % Equatorial elliptical, e_j < 0, r . v >= 0
            [Cart2Kep_UnitTest.defaultR, sqrt(2) / 2, 0, 0, 135, 225] ... % Equatorial elliptical, e_j >= 0, r . v < 0
            [Cart2Kep_UnitTest.defaultR, sqrt(2) / 2, 0, 0, 225, 225] ... % Equatorial elliptical, e_j < 0, r . v < 0
        };
            
    end
    
    methods (Test)
        % Tests invalid Cartesian elements (not a 1x6 or 6x1 vector)
        function testInvalidCartesianElements (testCase)
            testCase.verifyError(@() Cart2Kep([], 'True', 'Deg'), 'Cart2Kep:InvalidCartesianElements'); % Too few Cartesian elements
            testCase.verifyError(@() Cart2Kep([1, 1, 1, 1, 1, 1, 1], 'True', 'Deg'), 'Cart2Kep:InvalidCartesianElements'); % Too many Cartesian elements
            testCase.verifyError(@() Cart2Kep(['x', 'y', 'z', 'vx', 'vy', 'vz'], 'True', 'Deg'), 'Cart2Kep:InvalidCartesianElements'); % Invalid Cartesian elements
            testCase.verifyError(@() Cart2Kep([1, 1, 1; 1, 1, 1], 'True', 'Deg'), 'Cart2Kep:InvalidCartesianElements'); % Incorrect format
        end
        
        % Tests zero velocity vector
        function testZeroVelocity (testCase)
            testCase.verifyError(@() Cart2Kep([Cart2Kep_UnitTest.defaultR, Cart2Kep_UnitTest.defaultR, Cart2Kep_UnitTest.defaultR, 0, 0, 0], 'True', 'Deg'), 'Cart2Kep:ZeroVelocity');
        end
        
        % Tests invalid position (inside Earth)
        function testInvalidPosition (testCase)
            testCase.verifyError(@() Cart2Kep([0, 0, 0, 5, 5, 5], 'True', 'Deg'), 'Cart2Kep:InvalidPosition')
        end
        
        % Tests invalid anomaly types (not 'true,' 'mean,' or 'eccentric'
        function testInvalidAnomType (testCase, ver)
            testCase.verifyError(@() Cart2Kep(ver, 'False', 'Deg'), 'Cart2Kep:InvalidAnomType');
        end
        
       % Tests invalid unit types (not 'deg' or 'rad')
       function testInvalidUnitTypes (testCase, ver)
           testCase.verifyError(@() Cart2Kep(ver, 'True', 'Grad'), 'Cart2Kep:InvalidUnitType');
       end
       
       % Tests parabolic orbit (e = 1)
       function testParabolicOrbit (testCase, parabolic)
           testCase.verifyError(@() Cart2Kep(parabolic, 'True', 'Deg'), 'Cart2Kep:ParabolicOrbit');
       end
       
       % Tests parabolic orbit (e > 1)
       function testHyperbolicOrbit (testCase, hyperbolic)
           testCase.verifyError(@() Cart2Kep(hyperbolic, 'True', 'Deg'), 'Cart2Kep:HyperbolicOrbit');
       end
       
       % Tests rectilinear orbit (r x v = 0)
       function testRectilinearOrbit (testCase, rectilinear)
           testCase.verifyError(@() Cart2Kep(rectilinear, 'True', 'Deg'), 'Cart2Kep:RectilinearOrbit');
       end
       
       % Tests warnings for special case orbits
       function testSpecialCaseOrbitWarnings (testCase, equatorialCircular, inclinedCircular, equatorialElliptical)
           warning('on');
           testCase.verifyWarning(@() Cart2Kep(equatorialCircular, 'True', 'Deg'), 'Cart2Kep:EquatorialCircularOrbit'); % Equatorial circular orbit
           testCase.verifyWarning(@() Cart2Kep(inclinedCircular, 'True', 'Deg'), 'Cart2Kep:InclinedCircularOrbit'); % Inclined circular orbit
           testCase.verifyWarning(@() Cart2Kep(equatorialElliptical, 'True', 'Deg'), 'Cart2Kep:GeneralEquatorialOrbit'); % Equatorial elliptical orbit
       end
       
        % Tests Keplerian > Cartesian > Keplerian > Cartesian > Keplerian
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
                out = Cart2Kep(Kep2Cart(Cart2Kep(Kep2Cart(vec), 'Mean', 'Rad')), 'Mean', 'Rad');
                
                testCase.verifyEqual(out, vec, 'RelTol', 1E-5, 'AbsTol', 1E-5);
            end
        end
    end
    
    methods (Test, ParameterCombination = 'sequential')
        % Tests zeroes in Cartesian elements input
        function testCartesianElementZeroes (testCase, zeroVers, zeroResults)
            testCase.verifyEqual(Cart2Kep(zeroVers, 'True', 'Deg'), zeroResults, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end

        % Tests special case orbits
        function testSpecialCaseOrbits(testCase, specialVers, specialResults)
            testCase.verifyEqual(Cart2Kep(specialVers, 'True', 'Deg'), specialResults, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests against known case in degrees
        function testDegCases (testCase, degResVec, anomTypes)
            testCase.verifyEqual(Cart2Kep([6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341], anomTypes, 'Deg'), degResVec, 'RelTol', 1E-4, 'AbsTol', 1E-4);
        end

        % Tests against known cases in radians
        function testRadCases (testCase, radResVec, anomTypes)
            testCase.verifyEqual(Cart2Kep([6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341], anomTypes, 'Rad'), radResVec, 'RelTol', 1E-4, 'AbsTol', 1E-4);
        end
        
        % Tests Keplerian > Cartesian > Keplerian for special cases
        function testSpecialTwoWay (testCase, specialVers)
            testCase.verifyEqual(Kep2Cart(Cart2Kep(specialVers, 'Mean', 'Rad')), specialVers, 'RelTol', 1E-5, 'AbsTol', 1E-5);
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
% E. White       | 06-30-2023 | Added test for invalid input types and for
%                               conversion to Keplerian and back, added
%                               remaining edge case tests, added more
%                               two-way unit tests, made tests make use of
%                               default radius variables

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        convert_equinoctial_to_cartesian_UnitTest < matlab.unittest.TestCase
% convert_equinoctial_to_cartesian_UnitTest - Unit test for convert_equinoctial_to_cartesian
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
        maxR = 43000;
        minR = 6400;
    end
    
    properties (TestParameter)
        ver = {[0.000618165872034742, -0.186939977644424, 0.0251358447064123, -0.206130724502682, -2.36641891694603, -1.08358461740347, 0]};
        verRes = {[6513.711946044443, 6882.748803674651, 6438.350633995484, 4.902960531245754, -1.518827797965108, -1.983803852683226]};
    end
    
    methods (Test)
        % Tests default parameters
        function testDefaultParameters (testCase, ver)
            [actRvec,actVvec,actX,actY,actXdot,actYdot,actF,actCF,actSF] = convert_equinoctial_to_cartesian(ver(1), ver(2), ver(3), ver(4), ver(5), ver(6), ver(7));
            [expRvec,expVvec,expX,expY,expXdot,expYdot,expF,expCF,expSF] = convert_equinoctial_to_cartesian(ver(1), ver(2), ver(3), ver(4), ver(5), ver(6), ver(7), 1, 3.986004418e5, 100*eps(2*pi), 100, 2);
            act = [actRvec;actVvec;actX;actY;actXdot;actYdot;actF;actCF;actSF];
            exp = [expRvec;expVvec;expX;expY;expXdot;expYdot;expF;expCF;expSF];
            testCase.verifyEqual(act, exp);
        end
        
        % Tests against known verification case
        function testKnownCase(testCase, ver, verRes)
            [r, v] = convert_equinoctial_to_cartesian(ver(1), ver(2), ver(3), ver(4), ver(5), ver(6), ver(7));
            testCase.verifyEqual([r; v]', verRes, 'RelTol', 1E-5, 'AbsTol', 1E-5);
        end
        
        % Tests Equinoctial > Cartesian > Equinoctial > Cartesian >
        % Equinoctial
        function testTwoWay (testCase)
            for i = 1:1000
                r_1 = (convert_equinoctial_to_cartesian_UnitTest.maxR - convert_equinoctial_to_cartesian_UnitTest.minR) * rand + convert_equinoctial_to_cartesian_UnitTest.minR;
                r_2 = (convert_equinoctial_to_cartesian_UnitTest.maxR - convert_equinoctial_to_cartesian_UnitTest.minR) * rand + convert_equinoctial_to_cartesian_UnitTest.minR;

                a = (r_1 + r_2) / 2;
                e = rand;

                r_p = a * (1 - e);

                while r_p < convert_equinoctial_to_cartesian_UnitTest.minR
                    e = rand;
                    r_p = a * (1 - e);
                end

                vec = [a, e, pi * rand, 2 * pi * rand, 2 * pi * rand, 2 * pi * rand];
                cart = Kep2Cart(vec);
                [~,n,af,ag,chi,psi,lM] = convert_cartesian_to_equinoctial(cart(1:3), cart(4:6));
                eq = [n,af,ag,chi,psi,lM];
                [r, v] = convert_equinoctial_to_cartesian(n,af,ag,chi,psi,lM, 0);
                [~,n,af,ag,chi,psi,lM] = convert_cartesian_to_equinoctial(r, v);
                [r, v] = convert_equinoctial_to_cartesian(n,af,ag,chi,psi,lM, 0);
                [~,n,af,ag,chi,psi,lM] = convert_cartesian_to_equinoctial(r, v);
                out = [n,af,ag,chi,psi,lM];

                testCase.verifyEqual(eq, out, 'RelTol', 1E-5, 'AbsTol', 1E-5);
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

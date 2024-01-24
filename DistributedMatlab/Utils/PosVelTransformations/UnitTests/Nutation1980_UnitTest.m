classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        Nutation1980_UnitTest < matlab.unittest.TestCase
% Nutation1980_UnitTest - Unit test for Nutation1980
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Jul 2023;  Latest update: Jul 2023
%
% ----------------- BEGIN CODE -----------------
    
    properties (TestParameter)
        nutation4Terms = {[ ...
            0, 0, 0,  0, 1, -171996, -174.2, 92025,  8.9; ...
            0, 0, 2, -2, 2,  -13187,   -1.6,  5736, -3.1; ...
            0, 0, 2,  0, 2,   -2274,   -0.2,   977, -0.5; ...
            0, 0, 0,  0, 2,    2062,    0.2,  -895,  0.5 ...
            ]};
    end

    methods (Test)
        % Tests invalid nutation model
        function testInvalidNutationModel (testCase)
            testCase.verifyError(@() Nutation1980('678terms'), 'Nutation1980:InvalidModel');
        end
        
        % Tests first four nutation terms against Vallado (2004)
        function testFirstFourTerms (testCase, nutation4Terms)
            [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] = Nutation1980('4terms');
            
            act = [Li,Lpi,Fi,Di,OMi,S1i * 1E4,S2i * 1E4,C1i * 1E4,C2i * 1E4];
            exp = nutation4Terms;
            
            testCase.verifyEqual(act, exp, 'RelTol', 1E-10, 'AbsTol', 1E-10);
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
% E. White       | 07-12-2023 | Initial development

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

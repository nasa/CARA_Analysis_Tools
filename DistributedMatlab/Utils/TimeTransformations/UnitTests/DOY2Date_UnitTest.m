classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        DOY2Date_UnitTest < matlab.unittest.TestCase
% DOY2DATE_UnitTest - Unit test for DOY2Date
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
    
    properties (TestParameter)
        leapYears = {[2000; 2004], 2016};
        nonLeapYears = {[2003; 2100], 2023};
    end
    
    methods (Test)
        % Tests empty inputs
        function testEmptyInputs (testCase)
            testCase.verifyError(@() DOY2Date([], []), 'DOY2Date:InputMustBeNonempty');
        end
        
        % Tests column vector check
        function testColumnVectors (testCase)
            testCase.verifyError(@() DOY2Date([100, 200], [2003, 2010]), 'DOY2Date:InputsMustBeColumnVectors');
        end
        
        % Tests inputs of different sizes
        function testInputMismatch (testCase)
            testCase.verifyError(@() DOY2Date([100; 200; 300], [2003; 2010]), 'DOY2Date:InputSizesMustMatch');
        end
        
        % Tests non-numerical inputs for day and year
        function testNonNumericalInputs (testCase)
            testCase.verifyError(@() DOY2Date('a', 'b'), 'DOY2Date:InputsMustBeNumeric');
        end
        
        % Tests non-positive inputs for day and year
        function testNonPositiveInputs (testCase)
            testCase.verifyError(@() DOY2Date(0, 0), 'DOY2Date:InputsMustBePositive');
        end
        
        % Tests DOYs greater than 366
        function testDOYLimit (testCase)
            testCase.verifyError(@() DOY2Date(400, 2000), 'DOY2Date:DOYLimitExceeded');
        end
        
        % Tests for leap days used during leap years
        %
        % This test checks for both "normal" leap years (every four years)
        % and "special" leap years (every 400, since other multiples of 100
        % are not leap years)
        function testLeapYearInput (~, leapYears)
            DOY2Date(366 .* ones(size(leapYears)), leapYears);
        end
        
        % Tests for leap days used during non-leap years
        %
        % This test checks for both "normal" non-leap years
        % and "special" non-leap years (every 100, except for multiples of
        % 400)
        function testNonLeapYearInput (testCase, nonLeapYears)
            testCase.verifyError(@() DOY2Date(366 .* ones(size(nonLeapYears)), nonLeapYears), 'DOY2Date:InputYearNotLeap');
        end
    end
    
    methods (Test, ParameterCombination = 'sequential') 
        % Tests day sixty during leap years (February 29th)
        function testLeapSixty (testCase, leapYears)
            for i = 1:length(leapYears)
                [~, actResult] = DOY2Date(60, leapYears(i));
                expResult = [leapYears(i), 2, 29, 0, 0, 0];
                testCase.verifyEqual(actResult, expResult);
            end
        end
        
        % Tests day sixty during non-leap years (March 1st)
        function testNonLeapSixty (testCase, nonLeapYears)
            for i = 1:length(nonLeapYears)
                [~, actResult] = DOY2Date(60, nonLeapYears(i));
                expResult = [nonLeapYears(i), 3, 1, 0, 0, 0];
                testCase.verifyEqual(actResult, expResult);
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
% E. White       | 06-07-2023 |  Initial development

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        getCcdsTimeFormat_UnitTest < matlab.unittest.TestCase
% getCcsdsTimeFormat_UnitTest - Unit test for getCcsdsTimeFormat
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
    
    methods (Test)
        % Tests for missing "T" separator
        function testMissingTSeparator (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-06-11:12:55'), 'getCcsdsTimeFormat:NoTSeparatorFound');
        end
        
        % Tests for excess "T" separators
        function testExcessTSeparators (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-06T11:16:43T'), 'getCcsdsTimeFormat:MoreThanOneTSeparatorFound');
        end
        
        % Tests for invalid date formats
        function testInvalidDateFormat (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('06/06/23T11:18:01'), 'getCcsdsTimeFormat:InvalidDateFormat');
        end
        
        % Tests for invalid years
        function testInavlidYears (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('yyyy-06-06T04:18:47'), 'getCcsdsTimeFormat:InvalidDateFormat');
        end
        
        % Tests for invalid days and months
        function testInvalidDaysAndMonths (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-mm-06T04:19:29'), 'getCcsdsTimeFormat:InvalidDateFormat'); % Invalid month
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-ddT04:20:39'), 'getCcsdsTimeFormat:InvalidDateFormat'); % Invalid day
        end
        
        % Test for invalid hours, minutes, and seconds
        function testInvalidHMS (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-06THH:33:50'), 'getCcsdsTimeFormat:InvalidDateFormat'); % Invalid hour
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-06T04:MM:11'), 'getCcsdsTimeFormat:InvalidDateFormat'); % Invalid minute
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-06T04:34:SS'), 'getCcsdsTimeFormat:InvalidDateFormat'); % Invalid second
        end
        
        % Tests for invalid fractional second
        function testInvalidFraction (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-06T04:27:35.1Zs'), 'getCcsdsTimeFormat:InvalidDateFormat');
        end
        
        % Tests for excess second fraction separators
        function testExcessFractionSeparators (testCase)
            testCase.verifyWarning(@() getCcsdsTimeFormat('2023-06-06T11:19:41..53'), 'getCcsdsTimeFormat:MoreThanOneFractionSeparatorFound');
        end
        
        % Tests for decimal place with no trailing figures
        function testEmptyFraction (testCase)
            testCase.verifyTrue(endsWith(getCcsdsTimeFormat('2023-06-06T02:11:19.'), 'S'));
        end
        
        % Tests yyyy-DDDTHH:MM:SS.FZ format
        function testyyyyDDD (testCase)
            actResult = getCcsdsTimeFormat('2023-158T08:35:51.5Z');
            expResult = 'yyyy-DDDTHH:MM:SS.FZ';
            testCase.verifyEqual(actResult, expResult);
        end
        
        % Tests yyyy-mm-ddTHH:MM:SS.FZ format
        function testyyyymmdd (testCase)
            actResult = getCcsdsTimeFormat('2023-06-07T08:38:19.3Z');
            expResult = 'yyyy-mm-ddTHH:MM:SS.FZ';
            testCase.verifyEqual(actResult, expResult);
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

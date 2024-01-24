classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..'), ...
         matlab.unittest.fixtures.PathFixture('../../LoggingAndStringReporting')}) ...
        round_timestring_UnitTest < matlab.unittest.TestCase
% round_timestring_UnitTest - Unit test for round_timestring
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

    % Warning: some of the tests below may fail in future version of MATLAB
    % if floating point handling continues to change! In such a case, it
    % is recommended that you test the failing numbers in the console using
    % fprintf to high precision [i.e., fprintf('%.21f\n', <num>)] in both
    % the failing version of MATLAB and the latest tested working version.
    %
    % Latest tested working version: MATLAB R2019b
    
    methods (Test)
        % Tests invalid date/time input
        function testInvalidInput (testCase)
            testCase.verifyError(@() round_timestring('2023-158T08:54:35.1Z'), 'round_timestring:InvalidInputFormat')            
        end
        
        % Tests invalid time format
        function testInvalidTimeFormat(testCase)
            testCase.verifyError(@() round_timestring('2023-06-07 25:33:59'), 'round_timestring:InvalidTimeFormat'); % Invalid hour (> 24)
            testCase.verifyError(@() round_timestring('2023-06-07 09:60:30'), 'round_timestring:InvalidTimeFormat'); % Invalid minute (>= 60)
            testCase.verifyError(@() round_timestring('2023-06-07 09:31:60'), 'round_timestring:InvalidTimeFormat'); % Invalid second (>= 60)
            testCase.verifyError(@() round_timestring('2023-06-07 -1:34:40'), 'round_timestring:InvalidTimeFormat'); % Invalid hour (< 0)
            testCase.verifyError(@() round_timestring('2023-06-07 09:-1:56'), 'round_timestring:InvalidTimeFormat'); % Invalid minute (< 0)
            testCase.verifyError(@() round_timestring('2023-06-07 09:35:-1'), 'round_timestring:InvalidTimeFormat'); % Invalid second (< 0)
            testCase.verifyError(@() round_timestring('2023-06-07 08:54:35.1Z'), 'round_timestring:InvalidTimeFormat'); % Invalid format (has time zone indicator)
        end
        
        % Tests invalid time unit
        function testInvalidTimeUnit (testCase)
            testCase.verifyError(@() round_timestring('2023-06-07 09:39:40', 'ds'), 'round_timestring:InvalidTimeUnit');
        end
        
        % Tests default rounding case (ms)
        function testDefaultTimeUnit (testCase)
            testCase.verifyEqual(round_timestring('2023-06-07 09:45:11.12345'), '2023-06-07 09:45:11.123');
        end
        
        % Tests cases where time units must round up
        function testWrapTimeUnits (testCase)
            testCase.verifyEqual(round_timestring('2023-06-07 12:02:11.9999', 'ms'), '2023-06-07 12:02:12.000'); % Seconds
            testCase.verifyEqual(round_timestring('2023-06-07 12:21:59.9', 's'), '2023-06-07 12:22:00'); % Minutes
            testCase.verifyEqual(round_timestring('2023-06-07 12:59:59', 'm'), '2023-06-07 13:00'); % Hours
            testCase.verifyEqual(round_timestring('2023-06-07 23:59:59', 'h'), '2023-06-08 00:00'); % Days
        end
        
        % Tests all rounding cases
        function testRoundUnits (testCase)
            testCase.verifyEqual(round_timestring('2023-06-07 12:29:33.1234', 'ms'), '2023-06-07 12:29:33.123'); % Milliseconds down
            testCase.verifyEqual(round_timestring('2023-06-07 12:30:17.1236', 'ms'), '2023-06-07 12:30:17.124'); % Milliseconds up
            testCase.verifyEqual(round_timestring('2023-06-07 12:30:55.25', 's'), '2023-06-07 12:30:55'); % Seconds down
            testCase.verifyEqual(round_timestring('2023-06-07 12:31:23.75', 's'), '2023-06-07 12:31:24'); % Seconds up
            testCase.verifyEqual(round_timestring('2023-06-07 12:32:06', 'm'), '2023-06-07 12:32'); % Minutes down
            testCase.verifyEqual(round_timestring('2023-06-07 12:32:37', 'm'), '2023-06-07 12:33'); % Minutes up
            testCase.verifyEqual(round_timestring('2023-06-07 12:11:55', 'h'), '2023-06-07 12:00'); % Hours down
            testCase.verifyEqual(round_timestring('2023-06-07 12:33:11', 'h'), '2023-06-07 13:00'); % Hours up
            testCase.verifyEqual(round_timestring('2023-06-07 05:11:55', 'd'), '2023-06-07'); % Days down
            testCase.verifyEqual(round_timestring('2023-06-07 12:34:32', 'd'), '2023-06-08'); % Days up
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
% E. White       | 06-08-2023 |  Fixed broken test due to floating-point
%                                precision errors beterrn R2016b and R2019b

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

% EventRate_UnitTest - Unit test for EventRate
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Dependencies:
%
% /params/test_mission.m
%
% =========================================================================
%
% Initial version: Dec 2024; Latest update: Oct 2025
%
% ----------------- BEGIN CODE -----------------
% Run the examples provided for function EventRate

classdef (SharedTestFixtures = {...
           matlab.unittest.fixtures.PathFixture('../'), ...
           matlab.unittest.fixtures.PathFixture('../src') ...
           matlab.unittest.fixtures.PathFixture('../params') ...
           matlab.unittest.fixtures.PathFixture('../../ProbabilityOfCollision') ...
           matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
           matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
           matlab.unittest.fixtures.PathFixture('../../Utils/General') ...
           matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations') ...
           matlab.unittest.fixtures.PathFixture('../../Utils/Plotting') ...
           }) EventRate_UnitTest < matlab.unittest.TestCase
    methods (Test)
        function test1(testCase)
            Accuracy            = 0.01; 
            exp_medianEventRate = 0.074;
            exp_95EventRate     = 0.4;
            exp_PcumNoRMMs      = 0.001278093353127;
            exp_PcumRMMs        = 0.001235388521227;
            % Calculate Scale Factors
            actual_results_struct  = EventRate('test_mission');
            actual_medianEventRate = actual_results_struct.MissionEventRate(1);
            actual_95EventRate     = actual_results_struct.MissionEventRate(3);
            actual_PcumNoRMMs      = actual_results_struct.CumulativePcNoRMMs(1);
            actual_PcumRMMs        = actual_results_struct.CumulativePcWithRMMs(1);

            % Verify Results
            testCase.verifyEqual(actual_medianEventRate, exp_medianEventRate,'RelTol',Accuracy);
            testCase.verifyEqual(actual_95EventRate,  exp_95EventRate,'RelTol',Accuracy);
            testCase.verifyEqual(actual_PcumNoRMMs, exp_PcumNoRMMs,'RelTol',Accuracy);
            testCase.verifyEqual(actual_PcumRMMs, exp_PcumRMMs,'RelTol',Accuracy);

            % Clean up output files
            delete(gcp)
            rmdir output s
            close all
        end
    end
end



% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% N. Ravago      | 2024-Dec-12 | Initial development.
% N. Ravago      | 2025-Oct-29 | Changed to use MATLAB unit test framework
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
classdef (SharedTestFixtures = { ...
    matlab.unittest.fixtures.PathFixture('../') ...
    matlab.unittest.fixtures.PathFixture('../src') ...
    matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
    matlab.unittest.fixtures.PathFixture('../../Utils/General') ...
    matlab.unittest.fixtures.PathFixture('../../Utils/LoggingAndStringReporting') ...
    matlab.unittest.fixtures.PathFixture('../../Utils/Plotting') ...
    }) ...
    EvaluateLightPollutionUnitTest < matlab.unittest.TestCase
	
% EvaluateLightPollutionUnitTest - Unit Test for EvaluateLightPollution
%
% Syntax: results = runtests('EvaluateLightPollutionUnitTest')
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------
	
    methods (Test)
        function testLightPollutionEvaluation(testCase)
            % Load expected outputs
            ExpectedOutputs = load('ExpectedOutputs.mat').ExpectedOutputs;
            
            % Run test case
            output = run_EvaluateLightPollution_Test();
            
            % Verify all outputs
            testCase.verifyEqual(output.EvalTable, ExpectedOutputs.EvalTable, 'AbsTol', 1e-10);
            testCase.verifyEqual(output.Fbright, ExpectedOutputs.Fbright, 'AbsTol', 1e-10);
            testCase.verifyEqual(output.LightPollutionLevel, ExpectedOutputs.LightPollutionLevel, 'AbsTol', 1e-10);
            testCase.verifyEqual(output.LightPollutionRisk, ExpectedOutputs.LightPollutionRisk);
            testCase.verifyEqual(output.RedesignRecommendedFlag, ExpectedOutputs.RedesignRecommendedFlag);
            testCase.verifyEqual(output.LowLatEvalTable, ExpectedOutputs.LowLatEvalTable, 'AbsTol', 1e-10);
            testCase.verifyEqual(output.MidLatEvalTable, ExpectedOutputs.MidLatEvalTable, 'AbsTol', 1e-10);
            testCase.verifyEqual(output.HighLatEvalTable, ExpectedOutputs.HighLatEvalTable, 'AbsTol', 1e-10);
            
            % Cleanup outputs
            cleanup
        end
    end
end

function Output = run_EvaluateLightPollution_Test
    % Set up test parameters
    params.New.Nc              = 1584; 
    params.New.Altitude_km     =  550; 
    params.New.Inclination_deg = 53.0; 
    params.New.Mzen50    = 6.12;
    params.New.Mzen05    = 6.67;
    params.New.Mzen95    = 4.68;
    params.New.MzenAltkm = 550;
    params.Evaluation.SDAPoints = [0 45];
    
    % Suppress plots and output files (not needed for unit test)
    params.Plotting.On = false;
    params.Output.WriteMetricVsSDATable = false; 
    params.Output.WriteEvalReport = false;       
    params.Output.WriteEvalTable = false;  
    
    % Execute EvaluateLightPollution
    Output = EvaluateLightPollution(params);
end

function cleanup
    clear ans example examples failures i ME paramPath printStatement;
    if exist('output', 'dir')
        rmdir('output','s');
    end
end

% ----------------- END OF CODE -------------------------------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ----------------------------------------
% Developer      |    Date     |     Description
% -------------------------------------------------------------------------
% J.Halpin       |   11-25-24  | Initial development of unit test
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Pc3D_Hall_Utils') ...
        matlab.unittest.fixtures.PathFixture('../SDMC_Utils') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/General') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/OrbitTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations')}) ...
        PcMultiStepWithSDMC_UnitTest < matlab.unittest.TestCase ...
% PcMultiStepWithSDMC_UnitTest - Unit test for PcMultiStep
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Oct 2025;  Latest update: Oct 2025
%
% ----------------- BEGIN CODE -----------------

    properties (TestParameter)
        dataDirectory = {'../../../DataFiles/PcTestCaseCDMs'};
        pcMethodDataFile = {'CARA_PcMethod_Test_Conjunctions.xlsx'};
    end
    
    methods (Test)

        function testPcMultiStepWithSDMCHighAccuracyTop10(testCase, dataDirectory, pcMethodDataFile)
            [conjID, numTrials, numWorkers, sdmcPc] = testCase.ReadDataFile(dataDirectory, pcMethodDataFile);
            for i = 1:10
                disp(['Running SDMC test case ' num2str(i) ' in high-accuracy mode']);
                % Run PcMultiStepWithSDMC in high-accuracy mode using the
                % same number of trials and workers used to generate the
                % reference data
                cdmFile = fullfile(dataDirectory,[conjID{i} '.cdm']);
                [r1, v1, C1, r2, v2, C2, HBR, params] = ConjDecoder(cdmFile);
                params.apply_covXcorr_corrections = false;
                params.ForceSDMCCalculation = true;
                params.SDMCParams.num_trials = numTrials(i);
                params.SDMCParams.num_workers = numWorkers(i);
                [~,out] = PcMultiStepWithSDMC(r1,v1,C1,r2,v2,C2,HBR,params);
                Pc = out.SDMCPc;
                
                % We should have an exact match in this run mode
                testCase.verifyEqual(Pc,sdmcPc(i));
            end
        end
        
        function testPcMultiStepWithDefaultAccuracy(testCase, dataDirectory, pcMethodDataFile)
            [conjID, numTrials, ~, sdmcPc] = testCase.ReadDataFile(dataDirectory, pcMethodDataFile);
            numConj = length(conjID);
            for i = 1:numConj
                disp(['Running SDMC test case ' num2str(i) ' in default-accuracy mode']);
                % Run PcMultiStepWithSDMC in default-accuracy mode, this
                % will be less accurate than the high-accuracy mode, but
                % will run much faster. We'll have to use a binomial prop
                % test to verify that the outputs are ok.
                cdmFile = fullfile(dataDirectory,[conjID{i} '.cdm']);
                [r1, v1, C1, r2, v2, C2, HBR, params] = ConjDecoder(cdmFile);
                params.apply_covXcorr_corrections = false;
                params.ForceSDMCCalculation = true;
                [~,out] = PcMultiStepWithSDMC(r1,v1,C1,r2,v2,C2,HBR,params);
                Pc = out.SDMCPc;
                localNumTrials = out.SDMCInfo.numTrials;
                
                % At a pValue > 1e-3 we have good evidence that the Monte
                % Carlo results represent the same distribution
                [~,pValue] = binomial_prop_test(Pc,localNumTrials,sdmcPc(i),numTrials(i));
                testCase.verifyGreaterThan(pValue,1e-3);
            end
        end
    end
    
    methods (Access = private)
        function [conjID, numTrials, numWorkers, sdmcPc] = ReadDataFile(testCase, dataDirectory, pcMethodDataFile) %#ok<INUSL>
            caraData = readtable(fullfile(dataDirectory,pcMethodDataFile));
            conjID = caraData.Conjunction_ID;
            numTrials = caraData.NtotSDMC;
            numWorkers = caraData.NumSDMCWorkers;
            sdmcPc = caraData.PcSDMC;
        end
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% L. Baars       | 10-06-2025 | Initial Development

% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

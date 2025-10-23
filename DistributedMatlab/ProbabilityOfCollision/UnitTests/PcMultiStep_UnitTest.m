classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Pc3D_Hall_Utils') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/General') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/OrbitTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations')}) ...
        PcMultiStep_UnitTest < matlab.unittest.TestCase ...
% PcMultiStep_UnitTest - Unit test for PcMultiStep
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
        sampleR1 = {[-9.841950433215101e+05 +3.932342044549424e+05 +6.991223682230414e+06]};
        sampleV1 = {[+4.883696742000000e+03 +5.689086045000000e+03 +3.665361590000000e+02]};
        sampleC1 = {[+4.976545641899520e+04 +5.787130862568278e+04 +3.370410320935015e+03 +1.137272273949272e+01 -4.325472616114674e+00 -8.009705480233521e+01; ...
                     +5.787130862568278e+04 +6.730377643610841e+04 +3.926542932121541e+03 +1.321992688238858e+01 -5.035560720747812e+00 -9.314985106902773e+01; ...
                     +3.370410320935015e+03 +3.926542932121541e+03 +2.461403197221289e+02 +7.586865834476763e-01 -3.077848629905763e-01 -5.434034460756914e+00; ...
                     +1.137272273949272e+01 +1.321992688238858e+01 +7.586865834476763e-01 +2.608186227148725e-03 -9.804181796720670e-04 -1.829751672999786e-02; ...
                     -4.325472616114674e+00 -5.035560720747812e+00 -3.077848629905763e-01 -9.804181796720670e-04 +3.895883508545853e-04 +6.968892326415779e-03; ...
                     -8.009705480233521e+01 -9.314985106902773e+01 -5.434034460756914e+00 -1.829751672999786e-02 +6.968892326415779e-03 +1.289253320300791e-01]};
        sampleR2 = {[-9.839696058965517e+05 +3.936845951174244e+05 +6.991219291625473e+06]};
        sampleV2 = {[+1.509562687000000e+03 +7.372938617000000e+03 -1.492509430000000e+02]};
        sampleC2 = {[+4.246862551076427e+04 +2.066374367781032e+05 -5.011108933888592e+03 +3.104606531932427e+01 -1.201093683199582e+01 -2.207975848324051e+02; ...
                     +2.066374367781032e+05 +1.005854717283451e+06 -2.434876491048039e+04 +1.510022508670080e+02 -5.850063541467530e+01 -1.074752763805685e+03; ...
                     -5.011108933888592e+03 -2.434876491048039e+04 +6.131274993037449e+02 -3.667147183233717e+00 +1.391769957262238e+00 +2.601457791444154e+01; ...
                     +3.104606531932427e+01 +1.510022508670080e+02 -3.667147183233717e+00 +2.272826228568773e-02 -8.778253314778023e-03 -1.613538091053610e-01; ...
                     -1.201093683199582e+01 -5.850063541467530e+01 +1.391769957262238e+00 -8.778253314778023e-03 +3.428801115804722e-03 +6.251148178133809e-02; ...
                     -2.207975848324051e+02 -1.074752763805685e+03 +2.601457791444154e+01 -1.613538091053610e-01 +6.251148178133809e-02 +1.148404222181769e+00]};
        sampleHBR = {20};
        
        dataDirectory = {'../../../DataFiles/PcTestCaseCDMs'};
        pcMethodDataFile = {'CARA_PcMethod_Test_Conjunctions.xlsx'};
    end
    
    methods (Test)
        function testPriDefaultCov(testCase, sampleR1, sampleV1, sampleR2, sampleV2, sampleC2, sampleHBR)
            import matlab.unittest.constraints.HasNaN;
            
            sigmaValSq = (6378135 * 10)^2;
            localC1 = [sigmaValSq 0 0 0 0 0;
                       0 sigmaValSq 0 0 0 0;
                       0 0 sigmaValSq 0 0 0;
                       zeros(3,6)];
            
            [actPc,actOut] = PcMultiStep(sampleR1,sampleV1,localC1,sampleR2,sampleV2,sampleC2,sampleHBR);
            
            testCase.verifyThat(actPc,HasNaN);
            testCase.verifySubstring(actOut.PcMethod,'Unable to calculate Pc due to invalid 3x3 covariance matrices:');
            testCase.verifySubstring(actOut.PcMethod,'Primary Cov >= Default Covariance Cutoff');
        end
        
        function testSecDefaultCov(testCase, sampleR1, sampleV1, sampleC1, sampleR2, sampleV2, sampleHBR)
            import matlab.unittest.constraints.HasNaN;
            
            sigmaValSq = (6378135 * 10)^2;
            localC2 = [sigmaValSq 0 0 0 0 0;
                       0 sigmaValSq 0 0 0 0;
                       0 0 sigmaValSq 0 0 0;
                       zeros(3,6)];
            
            [actPc,actOut] = PcMultiStep(sampleR1,sampleV1,sampleC1,sampleR2,sampleV2,localC2,sampleHBR);
            
            testCase.verifyThat(actPc,HasNaN);
            testCase.verifySubstring(actOut.PcMethod,'Unable to calculate Pc due to invalid 3x3 covariance matrices:');
            testCase.verifySubstring(actOut.PcMethod,'Secondary Cov >= Default Covariance Cutoff');
        end
        
        function testPriCovNegVariance(testCase, sampleR1, sampleV1, sampleR2, sampleV2, sampleC2, sampleHBR)
            import matlab.unittest.constraints.HasNaN;
            
            sigmaValSq = -1;
            localC1 = [sigmaValSq 0 0 0 0 0;
                       0 sigmaValSq 0 0 0 0;
                       0 0 sigmaValSq 0 0 0;
                       zeros(3,6)];
            
            [actPc,actOut] = PcMultiStep(sampleR1,sampleV1,localC1,sampleR2,sampleV2,sampleC2,sampleHBR);
            
            testCase.verifyThat(actPc,HasNaN);
            testCase.verifySubstring(actOut.PcMethod,'Unable to calculate Pc due to invalid 3x3 covariance matrices:');
            testCase.verifySubstring(actOut.PcMethod,'Primary Cov has negative sigma^2 value');
        end
        
        function testSecCovNegVariance(testCase, sampleR1, sampleV1, sampleC1, sampleR2, sampleV2, sampleHBR)
            import matlab.unittest.constraints.HasNaN;
            
            sigmaValSq = -1;
            localC2 = [sigmaValSq 0 0 0 0 0;
                       0 sigmaValSq 0 0 0 0;
                       0 0 sigmaValSq 0 0 0;
                       zeros(3,6)];
            
            [actPc,actOut] = PcMultiStep(sampleR1,sampleV1,sampleC1,sampleR2,sampleV2,localC2,sampleHBR);
            
            testCase.verifyThat(actPc,HasNaN);
            testCase.verifySubstring(actOut.PcMethod,'Unable to calculate Pc due to invalid 3x3 covariance matrices:');
            testCase.verifySubstring(actOut.PcMethod,'Secondary Cov has negative sigma^2 value');
        end
        
        function testPriCovZeroValues(testCase, sampleR1, sampleV1, sampleR2, sampleV2, sampleC2, sampleHBR)
            import matlab.unittest.constraints.HasNaN;
            
            sigmaValSq = 0;
            localC1 = [sigmaValSq 0 0 0 0 0;
                       0 sigmaValSq 0 0 0 0;
                       0 0 sigmaValSq 0 0 0;
                       zeros(3,6)];
            
            [actPc,actOut] = PcMultiStep(sampleR1,sampleV1,localC1,sampleR2,sampleV2,sampleC2,sampleHBR);
            
            testCase.verifyThat(actPc,HasNaN);
            testCase.verifySubstring(actOut.PcMethod,'Unable to calculate Pc due to invalid 3x3 covariance matrices:');
            testCase.verifySubstring(actOut.PcMethod,'Primary Cov contains position components with values of zero');
        end
        
        function testSecCovZeroValues(testCase, sampleR1, sampleV1, sampleC1, sampleR2, sampleV2, sampleHBR)
            import matlab.unittest.constraints.HasNaN;
            
            sigmaValSq = 0;
            localC2 = [sigmaValSq 0 0 0 0 0;
                       0 sigmaValSq 0 0 0 0;
                       0 0 sigmaValSq 0 0 0;
                       zeros(3,6)];
            
            [actPc,actOut] = PcMultiStep(sampleR1,sampleV1,sampleC1,sampleR2,sampleV2,localC2,sampleHBR);
            
            testCase.verifyThat(actPc,HasNaN);
            testCase.verifySubstring(actOut.PcMethod,'Unable to calculate Pc due to invalid 3x3 covariance matrices:');
            testCase.verifySubstring(actOut.PcMethod,'Secondary Cov contains position components with values of zero');
        end
        
        function testPcMultiStepDefaultPcCalc(testCase, dataDirectory, pcMethodDataFile)
            [conjID, Pc2D, Nc2D, Nc3D, ...
                ViolPc2D, ViolNc2D, ViolNc3D] = testCase.ReadDataFile(dataDirectory, pcMethodDataFile);
            numConj = length(conjID);
            for i = 1:numConj
                % Run PcMultiStep in the default mode where it will
                % determine the best Pc based on usage violation logic
                cdmFile = fullfile(dataDirectory,[conjID{i} '.cdm']);
                [r1, v1, C1, r2, v2, C2, HBR, params] = ConjDecoder(cdmFile);
                params.apply_covXcorr_corrections = false;
                [Pc,out] = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR,params);
                
                % Verify that PcMultiStep calculates Pc2D usage violations
                % the same way
                testCase.verifyEqual(out.AnyPc2DViolations,ViolPc2D(i) > 0);
                if ViolPc2D(i) > 0
                    % Only if Pc2D violations are found, verify that
                    % PcMultiStep calculates Nc2D usage violations the same
                    % way
                    testCase.verifyEqual(out.AnyNc2DViolations,ViolNc2D(i) > 0);
                    if ViolNc2D(i) > 0
                        % Only if Nc2D violations are found, verify that
                        % PcMultiSTep calcualtes Nc3D usage violations the
                        % same way
                        testCase.verifyEqual(out.AnyNc3DViolations,ViolNc3D(i) > 0);
                        % Regardless of whether or not Nc3D usage
                        % violations are found, the best that PcMultiStep
                        % can do is produce an Nc3D estimate, verify the
                        % values are the same
                        testCase.verifyEqual(Pc,Nc3D(i),'RelTol',1e-3);
                    else
                        % If there are no Nc2D usage violations, then
                        % PcMultiStep should have returned the Nc2D value
                        % as the Pc
                        testCase.verifyEqual(Pc,Nc2D(i),'RelTol',1e-3);
                    end
                else
                    % If there are no Pc2D usage violations, then
                    % PcMultiSTep should have returned the Pc2D value as
                    % the Pc
                    testCase.verifyEqual(Pc,Pc2D(i),'RelTol',1e-12);
                end
            end
        end
        
        function testPcMultiStepForcedPc2D(testCase, dataDirectory, pcMethodDataFile)
            [conjID, Pc2D] = testCase.ReadDataFile(dataDirectory, pcMethodDataFile);
            numConj = length(conjID);
            for i = 1:numConj
                % Run PcMultiStep to only allow Pc2D calculations
                cdmFile = fullfile(dataDirectory,[conjID{i} '.cdm']);
                [r1, v1, C1, r2, v2, C2, HBR, params] = ConjDecoder(cdmFile);
                params.apply_covXcorr_corrections = false;
                params.OnlyPc2DCalculation = true;
                Pc = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR,params);
                
                % Pc reported by PcMultiStep should always be Pc2D in this
                % mode
                testCase.verifyEqual(Pc,Pc2D(i),'RelTol',1e-12);
            end
        end
        
        function testPcMultiStepForcedNc2D(testCase, dataDirectory, pcMethodDataFile)
            [conjID, ~, Nc2D] = testCase.ReadDataFile(dataDirectory, pcMethodDataFile);
            numConj = length(conjID);
            for i = 1:numConj
                % Run PcMultiStep to force Nc2D calculations
                cdmFile = fullfile(dataDirectory,[conjID{i} '.cdm']);
                [r1, v1, C1, r2, v2, C2, HBR, params] = ConjDecoder(cdmFile);
                params.apply_covXcorr_corrections = false;
                params.ForceNc2DCalculation = true;
                params.PreventNc3DCalculation = true;
                Pc = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR,params);
                
                % Pc reported by PcMultiStep should always be Nc2D in this
                % mode
                testCase.verifyEqual(Pc,Nc2D(i),'RelTol',1e-3);
            end
        end
        
        function testPcMultiStepForcedNc3D(testCase, dataDirectory, pcMethodDataFile)
            [conjID, ~, ~, Nc3D] = testCase.ReadDataFile(dataDirectory, pcMethodDataFile);
            numConj = length(conjID);
            for i = 1:numConj
                % Run PcMultiStep to force Nc2D calculations
                cdmFile = fullfile(dataDirectory,[conjID{i} '.cdm']);
                [r1, v1, C1, r2, v2, C2, HBR, params] = ConjDecoder(cdmFile);
                params.apply_covXcorr_corrections = false;
                params.ForceNc3DCalculation = true;
                Pc = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR,params);
                
                % Pc reported by PcMultiStep should always be Nc2D in this
                % mode
                testCase.verifyEqual(Pc,Nc3D(i),'RelTol',1e-3);
            end
        end
    end
    
    methods (Access = private)
        function [conjID, Pc2D, Nc2D, Nc3D, ...
                ViolPc2D, ViolNc2D, ViolNc3D] = ReadDataFile(testCase, dataDirectory, pcMethodDataFile) %#ok<INUSL>
            caraData = readtable(fullfile(dataDirectory,pcMethodDataFile));
            conjID = caraData.Conjunction_ID;
            Pc2D = caraData.Pc2D;
            Nc2D = caraData.Nc2D;
            Nc3D = caraData.Nc3D;
            ViolPc2D = caraData.ViolationsPc2D;
            ViolNc2D = caraData.ViolationsNc2D;
            ViolNc3D = caraData.ViolationsNc3D;
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

classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/LoggingAndStringReporting')}) ...
        PcCircleWithConjData_UnitTest < matlab.unittest.TestCase
% PcCircleWithConjData_UnitTest - Unit test for PcCircleWithConjData
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
% 
% Dependencies:
%
%   Pc_TestSetup.m
%
% =========================================================================
%
% Initial version: Dec 2023;  Latest update: Dec 2023
%
% ----------------- BEGIN CODE -----------------

    properties (TestParameter)
        multiR1   = {[-6227006.55254931 -3000839.63874701 -193256.909906263]};
        multiV1   = {[ 2983.98386544377  -5973.1889388784  -3622.2628437168]};
        multiCov1 = {[ 23292.4814989832 -46240.9751251139 -28103.5247270571;
                      -46240.9751251139  91893.7305844327  55838.5828667566;
                      -28103.5247270571  55838.5828667566  33965.0600024593]};
        multiR2   = {[-6231009.26896273  -2993196.1481304 -176140.165174705]};
        multiV2   = {[-325.596820869216  994.651454187559  -7507.3405392103]};
        multiCov2 = {[ 342571508.869164 -617504480.346085  4596554106.60093 -617504480.346085  390773572.222554 -2821108281.48108  4596554106.60093 -2821108281.48108  20336098305.5778;
                       36891498.9314705 -120066974.82088   900253524.718121  -120066974.82088  390773572.222554 -2929990982.07445  900253524.718121 -2929990982.07445  21968889281.9715]};
        multiHBR  = {10};
        multiPc   = {[0; 2.51493954968373e-15]};
        outIsPosDef = {[true; true]};
        outIsRemediated = {[true; false]};
        outSemiMajorAxis = {[114965.623083965; 123210.123996274]};
        outSemiMinorAxis = {[0.001; 18.6451744964948]};
        outClockAngle = {[11.4506649906806; -0.387821480125896]};
        outHBR = {[10; 10]};
        outMissDistance = {[19168.4020642769; 19168.4020642769]};
        outx1Sigma = {[0.00503717210711541; 2753.92453514032]};
        outRadialSigma = {[6897817.0037583; 1033485.46179984]};
        outInTrackSigma = {[38217015.3423524; 52748084.4916239]};
        outCrossTrackSigma = {[139862333.73159; 140047402.670176]};
        outConditionNumPrimary = {[7790.99912705967; 7790.99912705967]};
        outConditionNumSecondary = {[34369361.251137; 68183662.8273948]};
        outConditionNumCombined = {[210715.946879989; 64424477.1230422]};
        outConditionNumProjected = {[22.7597725912773; 43667656.4505498]};
        outRelativePhaseAngle = {[69.3767735325628; 69.3767735325628]};
    end
    
    methods (Test)
        % Tests extra output values, one at a time
        function testExtraOutputsSingly (testCase, multiR1, multiV1, multiCov1, multiR2, multiV2, multiCov2, multiHBR, ...
                multiPc, outIsPosDef, outIsRemediated, outSemiMajorAxis, outSemiMinorAxis, outClockAngle, outHBR, outMissDistance, ...
                outx1Sigma, outRadialSigma, outInTrackSigma, outCrossTrackSigma, outConditionNumPrimary, outConditionNumSecondary, ...
                outConditionNumCombined, outConditionNumProjected, outRelativePhaseAngle)
            
            numTests = length(multiPc);
            r1 = multiR1;
            v1 = multiV1;
            cov1 = multiCov1;
            r2 = multiR2;
            v2 = multiV2;
            HBR = multiHBR;
            relTol = 1e-9;
            
            for i = 1:numTests
                cov2 = [multiCov2(i,1:3); multiCov2(i,4:6); multiCov2(i,7:9)];
                [Pc, out] = PcCircleWithConjData(r1,v1,cov1,r2,v2,cov2,HBR);
                testCase.verifyEqual(Pc,multiPc(i),'RelTol',relTol);
                testCase.verifyEqual(out.IsPosDef,outIsPosDef(i));
                testCase.verifyEqual(out.IsRemediated,outIsRemediated(i));
                testCase.verifyEqual(out.SemiMajorAxis,outSemiMajorAxis(i),'RelTol',relTol);
                testCase.verifyEqual(out.SemiMinorAxis,outSemiMinorAxis(i),'RelTol',relTol);
                testCase.verifyEqual(out.ClockAngle,outClockAngle(i),'RelTol',relTol);
                testCase.verifyEqual(out.HBR,outHBR(i));
                testCase.verifyEqual(out.MissDistance,outMissDistance(i),'RelTol',relTol);
                testCase.verifyEqual(out.x1Sigma,outx1Sigma(i),'RelTol',relTol);
                testCase.verifyEqual(out.RadialSigma,outRadialSigma(i),'RelTol',relTol);
                testCase.verifyEqual(out.InTrackSigma,outInTrackSigma(i),'RelTol',relTol);
                testCase.verifyEqual(out.CrossTrackSigma,outCrossTrackSigma(i),'RelTol',relTol);
                testCase.verifyEqual(out.CondNumPrimary,outConditionNumPrimary(i),'RelTol',relTol);
                testCase.verifyEqual(out.CondNumSecondary,outConditionNumSecondary(i),'RelTol',relTol);
                testCase.verifyEqual(out.CondNumCombined,outConditionNumCombined(i),'RelTol',relTol);
                testCase.verifyEqual(out.CondNumProjected,outConditionNumProjected(i),'RelTol',relTol);
                testCase.verifyEqual(out.RelativePhaseAngle,outRelativePhaseAngle(i),'RelTol',relTol);
            end
        end
        
        % Tests extra output values, vectorized
        function testExtraOutputsVectorized (testCase, multiR1, multiV1, multiCov1, multiR2, multiV2, multiCov2, multiHBR, ...
                multiPc, outIsPosDef, outIsRemediated, outSemiMajorAxis, outSemiMinorAxis, outClockAngle, outHBR, outMissDistance, ...
                outx1Sigma, outRadialSigma, outInTrackSigma, outCrossTrackSigma, outConditionNumPrimary, outConditionNumSecondary, ...
                outConditionNumCombined, outConditionNumProjected, outRelativePhaseAngle)
            
            numTests = length(multiPc);
            r1 = repmat(multiR1,numTests,1);
            v1 = repmat(multiV1,numTests,1);
            cov1 = reshape(permute(multiCov1,[2 1]),1,9);
            cov1 = repmat(cov1,numTests,1);
            r2 = repmat(multiR2,numTests,1);
            v2 = repmat(multiV2,numTests,1);
            cov2 = multiCov2;
            HBR = repmat(multiHBR,numTests,1);
            relTol = 1e-9;
            
            [Pc, out] = PcCircleWithConjData(r1,v1,cov1,r2,v2,cov2,HBR);
            testCase.verifyEqual(Pc,multiPc,'RelTol',relTol);
            testCase.verifyEqual(out.IsPosDef,outIsPosDef);
            testCase.verifyEqual(out.IsRemediated,outIsRemediated);
            testCase.verifyEqual(out.SemiMajorAxis,outSemiMajorAxis,'RelTol',relTol);
            testCase.verifyEqual(out.SemiMinorAxis,outSemiMinorAxis,'RelTol',relTol);
            testCase.verifyEqual(out.ClockAngle,outClockAngle,'RelTol',relTol);
            testCase.verifyEqual(out.HBR,outHBR);
            testCase.verifyEqual(out.MissDistance,outMissDistance,'RelTol',relTol);
            testCase.verifyEqual(out.x1Sigma,outx1Sigma,'RelTol',relTol);
            testCase.verifyEqual(out.RadialSigma,outRadialSigma,'RelTol',relTol);
            testCase.verifyEqual(out.InTrackSigma,outInTrackSigma,'RelTol',relTol);
            testCase.verifyEqual(out.CrossTrackSigma,outCrossTrackSigma,'RelTol',relTol);
            testCase.verifyEqual(out.CondNumPrimary,outConditionNumPrimary,'RelTol',relTol);
            testCase.verifyEqual(out.CondNumSecondary,outConditionNumSecondary,'RelTol',relTol);
            testCase.verifyEqual(out.CondNumCombined,outConditionNumCombined,'RelTol',relTol);
            testCase.verifyEqual(out.CondNumProjected,outConditionNumProjected,'RelTol',relTol);
            testCase.verifyEqual(out.RelativePhaseAngle,outRelativePhaseAngle,'RelTol',relTol);
        end
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------------------------------
% Developer      |    Date    |     Description
%--------------------------------------------------------------------------
% L. Baars       | 12-27-2023 | Initial development

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
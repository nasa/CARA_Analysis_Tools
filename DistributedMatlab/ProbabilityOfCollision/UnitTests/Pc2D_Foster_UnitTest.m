classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/LoggingAndStringReporting')}) ...
        Pc2D_Foster_UnitTest < matlab.unittest.TestCase
% Pc2D_Foster_UnitTest - Unit test for Pc2D_Foster
%
% =========================================================================
%
% Copyright (c) 2019-2025 United States Government as represented by the
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
% Initial version: Dec 2019;  Latest update: Oct 2025
%
% ----------------- BEGIN CODE -----------------

    properties (TestParameter)
        AlfanoHBR = {15, 4, 15, 15, 10, 10, 10, 4, 6, 6, 4};
        AlfanoTLimit = {21600, 21600, 21600, 21600, 1419, 1419, 1419, 10135, 10800, 21600, 1420};
        AlfanoAccuracy = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};                  
        AlfanoExpSolution = { 1.46749549E-01, ...
                              6.22226700E-03, ...
                              1.00351176E-01, ...
                              4.93234060E-02, ...
                              4.44873860E-02, ...
                              4.33545500E-03, ...
                              1.58147000E-04, ...
                              3.69480080E-02, ...
                              2.90146291E-01, ...
                              2.90146291E-01, ...
                              2.67202600E-03 };

        AlfanoFileNum = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11'};
        
        r1 = {[378.39559 4305.721887 5752.767554]};
        v1 = {[2.360800244 5.580331936 -4.322349039]};
        cov1 = {[ 44.5757544811362   81.6751751052616  -67.8687662707124;
                  81.6751751052616  158.453402956163  -128.616921644857;
                 -67.8687662707124 -128.616921644858   105.490542562701]};
        
        r2 = {[374.5180598 4307.560983 5751.130418]};
        v2 = {[-5.388125081 -3.946827739 3.322820358]};
        cov2 = {[ 2.31067077720423  1.69905293875632 -1.4170164577661;
                  1.69905293875632  1.24957388457206 -1.04174164279599;
                 -1.4170164577661  -1.04174164279599  0.869260558223714]};
             
        OmitronHBR = {0.020};
    end

    methods (Test)        
        function testCircularHBR(testCase, r1, v1, cov1, r2, v2, cov2, OmitronHBR) % Original Test Omitron Test Case 01 (Written by Dragan Plakalovic)
            % Circular Area
            expSolution = 2.70601573490125e-05;
            Accuracy = 0.001; 
            Tol = 1e-09;
            HBRType = 'circle';
            
            [actSolution] = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,OmitronHBR,Tol,HBRType);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testSquareHBR(testCase, r1, v1, cov1, r2, v2, cov2, OmitronHBR) % Original Test Omitron Test Case 02 (Written by Dragan Plakalovic)
            % Square Area
            expSolution = 3.44534649703584e-05;
            Accuracy = 0.001; 
            Tol = 1e-09;
            HBRType = 'square';
            
            [actSolution] = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,OmitronHBR,Tol,HBRType);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testEqAreaSquareHBR(testCase, r1, v1, cov1, r2, v2, cov2, OmitronHBR) % Original Test Omitron Test Case 03 (Written by Dragan Plakalovic)
            % Square Equivalent Area
            expSolution = 2.70601573490111e-05;
            Accuracy = 0.001; 
            Tol = 1e-09;
            HBRType = 'squareEquArea';
            
            [actSolution] = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,OmitronHBR,Tol,HBRType);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testNonPositiveDefinite(testCase) % Non-PD test Case
            HBR = 20;
            Accuracy = 1E-10; 
            expSolution = 0;

            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/OmitronTestCase_Test07_NonPDCovariance.cdm');
                
            % Check for non- Positive Definite Error Throw
            [Pc2D,~,IsPosDef,IsRemediated]=Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            testCase.verifyEqual(Pc2D,expSolution,'AbsTol',Accuracy);
            testCase.verifyTrue(IsPosDef);
            testCase.verifyTrue(IsRemediated);
        end
        
        % Tests HBR type that is not 'circle', 'square', or 'squareEquArea'
        function testInvlidHBRType(testCase, r1, v1, cov1, r2, v2, cov2, OmitronHBR)
            Tol = 1e-09;
            HBRType = 'ellipse';
            
            testCase.verifyError(@() Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,OmitronHBR,Tol,HBRType), 'Pc2D_Foster:InvalidHBRType');
        end
    end
    
    methods (Test, ParameterCombination = 'sequential')
        function testAlfano(testCase, AlfanoHBR, AlfanoAccuracy, AlfanoExpSolution, AlfanoFileNum)
            % Read tables of A09 at-TCA conjunction data
            pars.case_list = str2double(AlfanoFileNum);
            pars.data_path = '../../../DataFiles/AlfanoInputData/Alfano_2009_Test_Cases';
            conj = GetAlfanoTestCases(pars);
            r1_ECI = conj.X1(1:3)';
            v1_ECI = conj.X1(4:6)';
            C1_ECI = conj.C1;
            r2_ECI = conj.X2(1:3)';
            v2_ECI = conj.X2(4:6)';
            C2_ECI = conj.C2;
            
            actSolution = Pc2D_Foster(r1_ECI,v1_ECI,C1_ECI,r2_ECI,v2_ECI,C2_ECI,AlfanoHBR,1e-8,'circle');
                
            testCase.verifyEqual(actSolution,AlfanoExpSolution,'RelTol',AlfanoAccuracy);
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
% T. Lechtenberg | 12-13-2019 | Initial Development
% L. Baars       | 10-03-2022 | Changed Alfano test case file reading
% E. White       | 08-07-2023 | Added compliant documentation, greatly
%                               improved readability and decreased length
%                               of unit tests
% E. White       | 08-09-2023 | Changed unit tests to use correct Alfano 
%                               data
% L. Baars       | 10-06-2025 | Updated path fixtures for more consistent
%                               passing runs.

% =========================================================================
%
% Copyright (c) 2019-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

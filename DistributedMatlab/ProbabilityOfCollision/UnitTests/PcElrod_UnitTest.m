classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/LoggingAndStringReporting')}) ...
        PcElrod_UnitTest < matlab.unittest.TestCase
% PcElrod_UnitTest - Unit test for PcElrod
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
% Initial version: Dec 2019;  Latest update: Aug 2023
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
    end

    methods (Test)        
        function testOmitronCase1(testCase) % Original Test Omitron Test Case 01 (Written by Dragan Plakalovic)
            expSolution = 2.70601573490125e-05;
            Accuracy = 0.001; 
            
            r1      = [378.39559 4305.721887 5752.767554];
            v1      = [2.360800244 5.580331936 -4.322349039];
            r2      = [374.5180598 4307.560983 5751.130418];
            v2      = [-5.388125081 -3.946827739 3.322820358];
            cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
                     81.6751751052616  158.453402956163  -128.616921644857;
                     -67.8687662707124 -128.616921644858 105.490542562701];
            cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
                     1.69905293875632  1.24957388457206  -1.04174164279599;
                     -1.4170164577661  -1.04174164279599 0.869260558223714];
            HBR     = 0.020;
            
            [actSolution] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testOmitronCase2(testCase) % Original Test Omitron Test Case 02 (Written by Dragan Plakalovic)
            expSolution = 0.180736305849476;
            Accuracy = 0.001; 
            
            r1      = [-3239.128337196251   2404.575152356222   5703.228541709001];
            v1      = [-3.745768373154199   5.012339015927846  -4.231864565717194];
            r2      = [-3239.138264917246   2404.568320465936   5703.235605231182];
            v2      = [6.110192790100711   -1.767321407894830   4.140369261741708];
            cov1    = [ 0.342072996423899  -0.412677096778269   0.371500417511149
                       -0.412677096778269   0.609905946319294  -0.540401385544286
                        0.371500417511149  -0.540401385544286   0.521238634755377].*1e-3;
            cov2    = [ 0.028351300975134  -0.008204103437377   0.019253747960155
                       -0.008204103437377   0.002404377774847  -0.005586512197914
                        0.019253747960155  -0.005586512197914   0.013289250260317];
            HBR     = 0.020;
            
            [actSolution] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testNonPositiveDefinite(testCase) % Original Test Omitron Test Case 03 (Written by Dragan Plakalovic)
            % Calculation with non-positive definite covariance
            expSolution = 0;
            Accuracy = 0.001;
            
            r1      = [-6711923.276204407   2115987.639388162   596903.033311516];
            v1      = [-3633.024009735033  -8343.876216289187   1529.269073307224];
            r2      = [-6703689.997789211   2066645.266249834   601175.621383896];
            v2      = [-499.951545210059   -8271.763735723972  -3673.287518932014];
            cov1    = [10651.06865424844    25162.51806981900  -4653.47179860534
                      25162.51806981901   61873.55078822979  -11370.63814221593
                     -4653.47179860534   -11370.63814221593   2108.71780971058];
            cov2    = [  38324418857.201    409236977944.646    182366534533.406
                      409236977944.646    4369926131386.673   1947351616381.380
                      182366534533.406    1947351616381.380   867790027954.990];
            HBR     = 52.84;
            
            [actSolution] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testNPDNoRemediation(testCase) % Original Test Omitron Test Case 04 (Written by Dragan Plakalovic)
            % Calculation with non-positive definite covariance that can't be remediated
            expSolution = NaN;
            Accuracy = 0.001;
            
            r1      = [-6711923.276204407   2115987.639388162   596903.033311516];
            v1      = [-3633.024009735033  -8343.876216289187   1529.269073307224];
            r2      = [-6703689.997789211   2066645.266249834   601175.621383896];
            v2      = [-499.951545210059   -8271.763735723972  -3673.287518932014];
            cov1    = [1.065106865424844e-3   25162.51806981900  -4653.47179860534
                      25162.51806981901     61873.55078822979  -11370.63814221593
                     -4653.47179860534     -11370.63814221593   2108.71780971058];
            cov2    = [  38324418857.201    409236977944.646    182366534533.406
                      409236977944.646    4369926131386.673   1947351616381.380
                      182366534533.406    1947351616381.380   867790027954.990];
            HBR     = 0.05284;
            
            [actSolution] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testMultiCalcNPD(testCase) % Original Test Omitron Test Case 05 (Written by Dragan Plakalovic)
            % Multiple calculations with one having an NPD matrix that can't be remediated
            expSolution = [0.180736305849476; NaN];
            Accuracy = 0.001;
            
            r1      = [-3239.128337196251   2404.575152356222   5703.228541709001;
                     -6711923.276204407   2115987.639388162   596903.033311516];
            v1      = [-3.745768373154199   5.012339015927846  -4.231864565717194;
                     -3633.024009735033  -8343.876216289187   1529.269073307224];
            r2      = [-3239.138264917246   2404.568320465936   5703.235605231182;
                     -6703689.997789211   2066645.266249834   601175.621383896];
            v2      = [6.110192790100711   -1.767321407894830   4.140369261741708;
                     -499.951545210059   -8271.763735723972  -3673.287518932014];
            cov1    = [ 0.342072996423899e-3 -0.412677096778269e-3  0.371500417511149e-3 ...
                     -0.412677096778269e-3  0.609905946319294e-3 -0.540401385544286e-3 ...
                      0.371500417511149e-3 -0.540401385544286e-3  0.521238634755377e-3;
                      1.065106865424844e-3  25162.51806981900    -4653.47179860534 ...
                      25162.51806981901     61873.55078822979    -11370.63814221593 ...
                     -4653.47179860534     -11370.63814221593     2108.71780971058];
            cov2    = [ 0.028351300975134  -0.008204103437377   0.019253747960155 ...
                     -0.008204103437377   0.002404377774847  -0.005586512197914 ...
                      0.019253747960155  -0.005586512197914   0.013289250260317;
                      38324418857.201     409236977944.646    182366534533.406 ...
                      409236977944.646    4369926131386.673   1947351616381.380 ...
                      182366534533.406    1947351616381.380   867790027954.990 ];
            HBR     = [0.020; 0.05284];
            
            [actSolution] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        function testNegativeHBR(testCase) % Original Test Omitron Test Case 06 (Written by Dragan Plakalovic)
            % HBR is negative
            expSolution = 0;
            Accuracy = 0.001;
            
            r1      = [378.39559 4305.721887 5752.767554];
            v1      = [2.360800244 5.580331936 -4.322349039];
            r2      = [374.5180598 4307.560983 5751.130418];
            v2      = [-5.388125081 -3.946827739 3.322820358];
            cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
                     81.6751751052616  158.453402956163  -128.616921644857;
                     -67.8687662707124 -128.616921644858 105.490542562701];
            cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
                     1.69905293875632  1.24957388457206  -1.04174164279599;
                     -1.4170164577661  -1.04174164279599 0.869260558223714];
            HBR     = -0.020;
            
            [actSolution] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end

        function testNonPositiveDefinite2(testCase) % Non-Pos Definite test Case
            HBR             = 20;
            Accuracy        = 1E-10; 
            expSolution     = 0;
            
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/OmitronTestCase_Test07_NonPDCovariance.cdm');
                
            % Check for non-Positive Definite Error Throw
            testCase.verifyWarning(@()PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,[],3),...
                    'RemediateCovariance2x2:stdNPD')
            [Pc,~,IsPosDef,IsRemediated] = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,[],3);
            testCase.verifyEqual(Pc,expSolution,'AbsTol',Accuracy);
            testCase.verifyTrue(IsPosDef);
            testCase.verifyTrue(IsRemediated);
        end
        
        % Returns a value very similar to Pc3D_Hall (see testSimilarTo2D in
        % Pc3D_Hall_UnitTest for more)
        function testSimilarTo3D(testCase)
            expSolution = 1.0281653e-02;
            Accuracy = 0.001; 
            
             r1   = [-9.841950433215101e+05 +3.932342044549424e+05 +6.991223682230414e+06];
             v1   = [+4.883696742000000e+03 +5.689086045000000e+03 +3.665361590000000e+02];
             cov1 = [+4.976545641899520e+04 +5.787130862568278e+04 +3.370410320935015e+03 +1.137272273949272e+01 -4.325472616114674e+00 -8.009705480233521e+01; ...
                     +5.787130862568278e+04 +6.730377643610841e+04 +3.926542932121541e+03 +1.321992688238858e+01 -5.035560720747812e+00 -9.314985106902773e+01; ...
                     +3.370410320935015e+03 +3.926542932121541e+03 +2.461403197221289e+02 +7.586865834476763e-01 -3.077848629905763e-01 -5.434034460756914e+00; ...
                     +1.137272273949272e+01 +1.321992688238858e+01 +7.586865834476763e-01 +2.608186227148725e-03 -9.804181796720670e-04 -1.829751672999786e-02; ...
                     -4.325472616114674e+00 -5.035560720747812e+00 -3.077848629905763e-01 -9.804181796720670e-04 +3.895883508545853e-04 +6.968892326415779e-03; ...
                     -8.009705480233521e+01 -9.314985106902773e+01 -5.434034460756914e+00 -1.829751672999786e-02 +6.968892326415779e-03 +1.289253320300791e-01];
             r2   = [-9.839696058965517e+05 +3.936845951174244e+05 +6.991219291625473e+06];
             v2   = [+1.509562687000000e+03 +7.372938617000000e+03 -1.492509430000000e+02];
             cov2 = [+4.246862551076427e+04 +2.066374367781032e+05 -5.011108933888592e+03 +3.104606531932427e+01 -1.201093683199582e+01 -2.207975848324051e+02; ...
                     +2.066374367781032e+05 +1.005854717283451e+06 -2.434876491048039e+04 +1.510022508670080e+02 -5.850063541467530e+01 -1.074752763805685e+03; ...
                     -5.011108933888592e+03 -2.434876491048039e+04 +6.131274993037449e+02 -3.667147183233717e+00 +1.391769957262238e+00 +2.601457791444154e+01; ...
                     +3.104606531932427e+01 +1.510022508670080e+02 -3.667147183233717e+00 +2.272826228568773e-02 -8.778253314778023e-03 -1.613538091053610e-01; ...
                     -1.201093683199582e+01 -5.850063541467530e+01 +1.391769957262238e+00 -8.778253314778023e-03 +3.428801115804722e-03 +6.251148178133809e-02; ...
                     -2.207975848324051e+02 -1.074752763805685e+03 +2.601457791444154e+01 -1.613538091053610e-01 +6.251148178133809e-02 +1.148404222181769e+00];
             HBR     = 20;
            
            [actSolution] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        % Tests default parameter values
        function testDefaultParameters (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase01.cdm');
            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,15);
            expSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,15,64,0);
            
            testCase.verifyEqual(actSolution,expSolution);
        end
        
        % Tests unequal number of primary and secondary positions
        function testUnequalPositionCount (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase01.cdm');
            
            r2_J2K = [r2_J2K; r2_J2K];
            
            testCase.verifyError(@() PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,15), 'PcElrod:UnequalPositionCount');
        end
        
        % Tests HBR values that are not [1x1] or [nx1]
        function testInvalidHBRDimensions (testCase)
            [r1_J2K_1, v1_J2K_1, C1_J2K_1, r2_J2K_1, v2_J2K_1, C2_J2K_1] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase01.cdm');
            [r1_J2K_2, v1_J2K_2, C1_J2K_2, r2_J2K_2, v2_J2K_2, C2_J2K_2] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase01.cdm');
            
            r1_J2K = [r1_J2K_1; r1_J2K_2];
            v1_J2K = [v1_J2K_1; v1_J2K_2];
            C1_J2K(:, :, 1) = C1_J2K_1;
            C1_J2K(:, :, 2) = C1_J2K_2;
            r2_J2K = [r2_J2K_1; r2_J2K_2];
            v2_J2K = [v2_J2K_1; v2_J2K_2];
            C2_J2K(:, :, 1) = C2_J2K_1;
            C2_J2K(:, :, 2) = C2_J2K_2;
            
            testCase.verifyError(@() PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,[15 4]), 'PcElrod:InvalidHBRDimensions'); % [1xn] instead of [nx1]
            testCase.verifyError(@() PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,[15 15; 4 4]), 'PcElrod:InvalidHBRDimensions'); % [nxn] instead of [nx1]
        end
        
        % Tests negative HBR value warning
        function testNegativeHBRWarning (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase01.cdm');
            
            testCase.verifyWarning(@() PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,-15,64,1), 'PcElrod:NegativeHBR');
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
            
            actSolution = PcElrod(r1_ECI,v1_ECI,C1_ECI,r2_ECI,v2_ECI,C2_ECI,AlfanoHBR);
                
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
% T. Lechtenberg | 03-19-2020 | Changed PcElrod call to excplicitly utilize
%                               ChebyshevOrder and WarningLevel parameters
% L. Baars       | 10-03-2022 | Changed Alfano test case file reading and
%                               custom error parameter
% E. White       | 08-07-2023 | Added compliant documentation, greatly
%                               improved readability and decreased length
%                               of unit tests, added additional unit tests
%                               from verification cases given in PcElrod.m,
%                               added tests for warning/error handling,
%                               added test for default parameters
% E. White       | 08-09-2023 | Changed unit tests to use correct Alfano 
%                               data

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

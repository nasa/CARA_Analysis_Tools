classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/LoggingAndStringReporting')}) ...
        PcCircle_UnitTest < matlab.unittest.TestCase
% PcCircle_UnitTest - Unit test for PcCircle
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
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
% Initial version: Aug 2023;  Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------
    
    properties (TestParameter)
        AlfanoHBR = {15, 4, 15, 15, 10, 10, 10, 4, 6, 6, 4, 4};
        AlfanoAccuracy = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};                  
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
                              2.67202600E-03, ...
                              NaN};
        AlfanoFileNum = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'};
        
        omitronR1      = {[378.39559 4305.721887 5752.767554]};
        omitronV1      = {[2.360800244 5.580331936 -4.322349039]};
        omitronR2      = {[374.5180598 4307.560983 5751.130418]};
        omitronV2      = {[-5.388125081 -3.946827739 3.322820358]};
        omitronCov1    = {[ 44.5757544811362   81.6751751052616  -67.8687662707124;
                            81.6751751052616  158.453402956163  -128.616921644857;
                           -67.8687662707124 -128.616921644858   105.490542562701]};
        omitronCov2    = {[ 2.31067077720423  1.69905293875632  -1.4170164577661;
                            1.69905293875632  1.24957388457206  -1.04174164279599;
                           -1.4170164577661  -1.04174164279599   0.869260558223714]};
        omitronHBR     = {0.020};
    end

    methods (Test)        
        function testOmitronCase1a (testCase, omitronR1, omitronV1, omitronR2, omitronV2, omitronCov1, omitronCov2, omitronHBR) % Default-method output
            expSolution = 2.706023476569787e-05;
          
            actSolution = PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR);
          
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1E-10);
        end
        
        function testOmitronCase1b (testCase, omitronR1, omitronV1, omitronR2, omitronV2, omitronCov1, omitronCov2, omitronHBR) % Highest-accuracy output
            expSolution = 2.706023476569786e-05;
            
            params.EstimationMode = 1;
            actSolution = PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR, params);
          
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1E-10);
        end
        
        function testOmitronCase2 (testCase) % Highest-accuracy output
            expSolution = 1.807363058494765e-01;
            
            r1      = [-3239.128337196251   2404.575152356222   5703.228541709001];
            v1      = [-3.745768373154199   5.012339015927846  -4.231864565717194];
            r2      = [-3239.138264917246   2404.568320465936   5703.235605231182];
            v2      = [6.110192790100711   -1.767321407894830   4.140369261741708];
            cov1    = [ 0.342072996423899 -0.412677096778269  0.371500417511149
                       -0.412677096778269  0.609905946319294 -0.540401385544286
                        0.371500417511149 -0.540401385544286  0.521238634755377].*1e-3;
            cov2    = [ 0.028351300975134 -0.008204103437377  0.019253747960155
                       -0.008204103437377  0.002404377774847 -0.005586512197914
                        0.019253747960155 -0.005586512197914  0.013289250260317];
            HBR     = 0.020;
            
            params.WarningLevel = 1;
            actSolution = PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params);
          
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1E-10);
        end
        
        function testNPDCovariance (testCase) % Example with non-positive definite covariance
            expSolution = 0;
            
            r1      = [-6711923.276204407   2115987.639388162   596903.033311516];
            v1      = [-3633.024009735033  -8343.876216289187   1529.269073307224];
            r2      = [-6703689.997789211   2066645.266249834   601175.621383896];
            v2      = [-499.951545210059   -8271.763735723972  -3673.287518932014];
            cov1    = [10651.06865424844  25162.51806981900  -4653.47179860534
                       25162.51806981901  61873.55078822979 -11370.63814221593
                       -4653.47179860534 -11370.63814221593   2108.71780971058];
            cov2    = [ 38324418857.201  409236977944.646  182366534533.406
                       409236977944.646 4369926131386.673 1947351616381.380
                       182366534533.406 1947351616381.380  867790027954.990];
            HBR     = 52.84;
            
            params.WarningLevel = 1;
            actSolution = PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params);
            
            testCase.verifyWarning(@() PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params), 'PcCircle:NPDCovariance');
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1E-10);
        end
        
        function testTwoElementNPDCovariance (testCase) % Two-element example calculation with one having an NPD matrix
            expSolution = [1.807363058494765e-01; 0];
            
            r1      = [   -3239.128337196251    2404.575152356222   5703.228541709001;
                       -6711923.276204407    2115987.639388162    596903.033311516];
            v1      = [   -3.745768373154199     5.012339015927846   -4.231864565717194;
                       -3633.024009735033    -8343.876216289187    1529.269073307224];
            r2      = [   -3239.138264917246    2404.568320465936   5703.235605231182;
                       -6703689.997789211    2066645.266249834    601175.621383896];
            v2      = [   6.110192790100711    -1.767321407894830     4.140369261741708;
                       -499.951545210059    -8271.763735723972    -3673.287518932014];
            cov1    = [    0.342072996423899e-3     -0.412677096778269e-3      0.371500417511149e-3 ...
                          -0.412677096778269e-3      0.609905946319294e-3     -0.540401385544286e-3 ...
                           0.371500417511149e-3     -0.540401385544286e-3      0.521238634755377e-3;
                       10651.06865424844         25162.51806981900         -4653.47179860534 ...
                       25162.51806981901         61873.55078822979        -11370.63814221593 ...
                       -4653.47179860534        -11370.63814221593          2108.71780971058];
            cov2    = [          0.028351300975134            -0.008204103437377             0.019253747960155 ...
                                -0.008204103437377             0.002404377774847            -0.005586512197914 ...
                                 0.019253747960155            -0.005586512197914             0.013289250260317;
                       38324418857.201              409236977944.646              182366534533.406 ...
                       409236977944.646            4369926131386.673             1947351616381.380 ...
                       182366534533.406            1947351616381.380              867790027954.990 ];
            HBR     = [0.020; 52.84];
            
            params.WarningLevel = 1;
            actSolution = PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params);
            
            testCase.verifyWarning(@() PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params), 'PcCircle:NPDCovariance');
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1E-10);
        end
        
        function testMultipleElements (testCase) % Multiple-element calculation with HBR values spanning 0 to Inf
            expSolution = [0; ...
                            2.514939701718358e-15; ...
                            1.232312816404174e-05; ...
                            6.343318773183237e-03; ...
                            1.917884657463699e-02; ...
                            6.389821634646267e-02; ...
                            1.900663422510684e-01; ...
                            1.000000000000000e+00];

            r1 = [-6.227006552549310e+06    -3.000839638747010e+06    -1.932569099062626e+05];
            v1 = [ 2.983983865443772e+03    -5.973188938878403e+03    -3.622262843716804e+03];
            c1 = [ 2.329248149898317e+04    -4.624097512511387e+04    -2.810352472705708e+04; ...
                  -4.624097512511387e+04     9.189373058443272e+04     5.583858286675660e+04; ...
                  -2.810352472705708e+04     5.583858286675660e+04     3.396506000245927e+04];   
            r2 = [-6.231009268962734e+06    -2.993196148130401e+06    -1.761401651747054e+05];
            v2 = [-3.255968208692156e+02     9.946514541875586e+02    -7.507340539210305e+03];
            c2 = [ 3.689149893147051e+07    -1.200669748208798e+08     9.002535247181211e+08; ...
                  -1.200669748208798e+08     3.907735722225543e+08    -2.929990982074453e+09; ...
                   9.002535247181211e+08    -2.929990982074453e+09     2.196888928197153e+10];
            HBR = [0 1e1 1e2 1e3 3e3 1e4 3e4 Inf]'; NHBR = numel(HBR);
            r1rep = repmat(r1,[NHBR 1]);
            v1rep = repmat(v1,[NHBR 1]);
            c1rep = repmat(c1,[1 1 NHBR]);
            r2rep = repmat(r2,[NHBR 1]);
            v2rep = repmat(v2,[NHBR 1]);
            c2rep = repmat(c2,[1 1 NHBR]);
            
            params.WarningLevel = 1;
            actSolution = PcCircle(r1rep,v1rep,c1rep,r2rep,v2rep,c2rep,HBR,params);
            
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1E-10);
        end
        
        % Tests both invalid (error-causing) and insufficient
        % (warning-causing, leads to innacuracies) EstimationModes
        function testInvalidAndInsufficentEstimationMode (testCase, omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR)
            % Tests non-positive EstimationMode that is neither 0 nor -1
            params.EstimationMode = -2;
            testCase.verifyError(@() PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR, params), 'PcCircle:InvalidEstimationMode');
            
            % Tests non-integer EstimationMode
            params.EstimationMode = 25.5;
            testCase.verifyError(@() PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR, params), 'PcCircle:InvalidEstimationMode');
            
            % Tests invalid number of arguments
            testCase.verifyError(@() PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2), 'PcCircle:InvalidNumberOfArguments');
            
            % Tests EstimationMode < 16
            params.EstimationMode = 11;
            params.WarningLevel = 1;
            testCase.verifyWarning(@() PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR, params), 'PcCircle:InsufficientEstimationMode');
        end
        
        % Tests unequal number of primary and secondary positions
        function testUnequalPositionCount (testCase, omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR)
            omitronR2 = [omitronR2; omitronR2];
            
            testCase.verifyError(@() PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR), 'PcCircle:UnequalPositionCount');
        end
        
        % Tests HBR values that are not [1x1] or [nx1]
        function testInvalidHBRDimensions (testCase, omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR)
            r1 = [omitronR1; omitronR1];
            v1 = [omitronV1; omitronV1];
            cov1 = reshape(omitronCov1, [1, 9]);
            cov1 = [cov1; cov1];
            r2 = [omitronR2; omitronR2];
            v2 = [omitronV2; omitronV2];
            cov2 = reshape(omitronCov2, [1, 9]);
            cov2 = [cov2; cov2];
            
            HBR = [omitronHBR, omitronHBR];
            testCase.verifyError(@() PcCircle(r1, v1, cov1, r2, v2, cov2, HBR), 'PcCircle:InvalidHBRDimensions'); % [1xn] instead of [nx1]
            
            HBR = [HBR; HBR];
            testCase.verifyError(@() PcCircle(r1, v1, cov1, r2, v2, cov2, HBR), 'PcCircle:InvalidHBRDimensions'); % [nxn] instead of [nx1]
        end
        
        % Tests negative HBR value warning
        function testNegativeHBR (testCase, omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, omitronHBR)
            params.EstimationMode = 64;
            params.WarningLevel = 1;
            testCase.verifyWarning(@() PcCircle(omitronR1, omitronV1, omitronCov1, omitronR2, omitronV2, omitronCov2, -omitronHBR, params), 'PcCircle:NegativeHBR');
        end
        
        function testZeroMissRelVelSingly (testCase) % Tests regular conjunctions with zero miss distance cases one at a time
            expSolution = [1.807363058494765e-01;
                4.24233551288858e-07;
                4.2421193442297e-07;
                NaN;
                1.56178901385447e-06;
                NaN;
                1.27770802292292e-08];
            
            r1   = [-3239.128337196251   2404.575152356222   5703.228541709001];
            v1   = [-3.745768373154199   5.012339015927846  -4.231864565717194];
            r2   = [-3239.138264917246   2404.568320465936   5703.235605231182];
            v2   = [6.110192790100711   -1.767321407894830   4.140369261741708];
            cov1 = [ 0.342072996423899 -0.412677096778269  0.371500417511149 -0.412677096778269  0.609905946319294 -0.540401385544286 0.371500417511149 -0.540401385544286  0.521238634755377].*1e-3;
            cov2 = [ 0.028351300975134 -0.008204103437377  0.019253747960155 -0.008204103437377  0.002404377774847 -0.005586512197914 0.019253747960155 -0.005586512197914  0.013289250260317];
            HBR  = 0.020;
            
            params.WarningLevel = 1;
            actSolution = PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params);
            testCase.verifyEqual(actSolution, expSolution(1), 'RelTol', 1E-10);
            
            % Load and run 6 of the Alfano test cases
            for i = 1:6
                fileName = sprintf('../../../DataFiles/SampleCDMs/AlfanoTestCase%02d.cdm',i);
                [r1, v1, C1, r2, v2, C2] = Pc_TestSetup(fileName);
                % Modify cases 2 and 4 for zero miss distance
                if i == 2 || i == 4
                    r2 = r1;
                % Modify cases 3 and 5 for zero relative velocity
                elseif i == 3 || i == 5
                    v2 = v1;
                end
                HBR = 0.010;
                
                % Verify that warnings get thrown for individual test cases
                if i == 2 || i == 4
                    testCase.verifyWarning(@() PcCircle(r1,v1,C1,r2,v2,C2,HBR,params), 'PcCircle:ZeroMissDistance');
                elseif i == 3 || i == 5
                    testCase.verifyWarning(@() PcCircle(r1,v1,C1,r2,v2,C2,HBR,params), 'PcCircle:ZeroRelativeVelocity');
                end
                
                % Verify that no warnings are thrown when WarningLevel
                % isn't set
                testCase.verifyWarningFree(@() PcCircle(r1,v1,C1,r2,v2,C2,HBR));
                
                % Verify the actual solution values
                actSolution = PcCircle(r1,v1,C1,r2,v2,C2,HBR,params);
                testCase.verifyEqual(actSolution, expSolution(i+1), 'RelTol', 1E-10);
            end
        end
        
        function testZeroMissRelVelVectorized (testCase) % Tests regular conjunctions with zero miss distance cases using vectorized calculations
            expSolution = [1.807363058494765e-01;
                4.24233551288858e-07;
                4.2421193442297e-07;
                NaN;
                1.56178901385447e-06;
                NaN;
                1.27770802292292e-08];
            
            numCases = 7;
            r1 = nan(numCases,3);
            v1 = nan(numCases,3);
            cov1 = nan(numCases,9);
            r2 = nan(numCases,3);
            v2 = nan(numCases,3);
            cov2 = nan(numCases,9);
            HBR = nan(numCases,1);
            
            r1(1,:)   = [-3239.128337196251   2404.575152356222   5703.228541709001];
            v1(1,:)   = [-3.745768373154199   5.012339015927846  -4.231864565717194];
            r2(1,:)   = [-3239.138264917246   2404.568320465936   5703.235605231182];
            v2(1,:)   = [6.110192790100711   -1.767321407894830   4.140369261741708];
            cov1(1,:) = [ 0.342072996423899 -0.412677096778269  0.371500417511149 -0.412677096778269  0.609905946319294 -0.540401385544286 0.371500417511149 -0.540401385544286  0.521238634755377].*1e-3;
            cov2(1,:) = [ 0.028351300975134 -0.008204103437377  0.019253747960155 -0.008204103437377  0.002404377774847 -0.005586512197914 0.019253747960155 -0.005586512197914  0.013289250260317];
            HBR(1)    = 0.020;
            
            % Load 6 of the Alfano test cases
            for i = 1:6
                fileName = sprintf('../../../DataFiles/SampleCDMs/AlfanoTestCase%02d.cdm',i);
                [r1a, v1a, C1a, r2a, v2a, C2a] = Pc_TestSetup(fileName);
                % Modify cases 2 and 4 for zero miss distance
                if i == 2 || i == 4
                    r2a = r1a;
                % Modify cases 3 and 5 for zero relative velocity
                elseif i == 3 || i == 5
                    v2a = v1a;
                end
                r1(i+1,:)   = r1a;
                v1(i+1,:)   = v1a;
                cov1(i+1,:) = [C1a(1,1:3) C1a(2,1:3) C1a(3,1:3)];
                r2(i+1,:)   = r2a;
                v2(i+1,:)   = v2a;
                cov2(i+1,:) = [C2a(1,1:3) C2a(2,1:3) C2a(3,1:3)];
                HBR(i+1)    = 0.010;
            end
            
            % Verify that warnings get thrown
            params.WarningLevel = 1;
            testCase.verifyWarning(@() PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params), 'PcCircle:ZeroMissDistance');
            testCase.verifyWarning(@() PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params), 'PcCircle:ZeroRelativeVelocity');
            
            % Verify that no warnings are thrown when WarningLevel isn't
            % set
            testCase.verifyWarningFree(@() PcCircle(r1,v1,cov1,r2,v2,cov2,HBR));
            
            % Verify the actual solution values
            actSolution = PcCircle(r1,v1,cov1,r2,v2,cov2,HBR);
            testCase.verifyEqual(actSolution, expSolution, 'RelTol', 1E-10);
        end

        % Tests all Alfano test cases in vectorized mode
        function testAlfanoVectorized(testCase)
            numRows = length(testCase.AlfanoHBR);
            HBR = cell2mat(testCase.AlfanoHBR)';
            accuracy = cell2mat(testCase.AlfanoAccuracy)';
            expSolution = cell2mat(testCase.AlfanoExpSolution)';
            fileNum = testCase.AlfanoFileNum;
            r1 = nan(numRows,3);
            v1 = nan(numRows,3);
            C1 = nan(numRows,9);
            r2 = nan(numRows,3);
            v2 = nan(numRows,3);
            C2 = nan(numRows,9);
            for i = 1:numRows
                params.case_list = str2double(fileNum{i});
                params.data_path = '../../../DataFiles/AlfanoInputData/Alfano_2009_Test_Cases';
                conj = GetAlfanoTestCases(params);
                r1(i,:) = conj.X1(1:3)';
                v1(i,:) = conj.X1(4:6)';
                C1(i,:) = [conj.C1(1,1) conj.C1(2,1) conj.C1(3,1) conj.C1(1,2) conj.C1(2,2) conj.C1(3,2) conj.C1(1,3) conj.C1(2,3) conj.C1(3,3)];
                r2(i,:) = conj.X2(1:3)';
                v2(i,:) = conj.X2(4:6)';
                C2(i,:) = [conj.C2(1,1) conj.C2(2,1) conj.C2(3,1) conj.C2(1,2) conj.C2(2,2) conj.C2(3,2) conj.C2(1,3) conj.C2(2,3) conj.C2(3,3)];
            end

            actSolution = PcCircle(r1,v1,C1,r2,v2,C2,HBR);
            for i = 1:numRows
                testCase.verifyEqual(actSolution(i),expSolution(i),'RelTol',accuracy(i));
            end
        end

        % Tests all Alfano test cases in vectorized mode, interspersing the
        % special cases (i.e., zero miss distance and zero relative
        % velocity) throughout the calculation
        function testAlfanoVectorizedWithSpecialCases(testCase)
            numRows = length(testCase.AlfanoHBR);
            HBR = cell2mat(testCase.AlfanoHBR)';
            accuracy = cell2mat(testCase.AlfanoAccuracy)';
            expSolution = cell2mat(testCase.AlfanoExpSolution)';
            fileNum = testCase.AlfanoFileNum;
            r1 = nan(numRows,3);
            v1 = nan(numRows,3);
            C1 = nan(numRows,9);
            r2 = nan(numRows,3);
            v2 = nan(numRows,3);
            C2 = nan(numRows,9);
            for i = 1:numRows
                params.case_list = str2double(fileNum{i});
                params.data_path = '../../../DataFiles/AlfanoInputData/Alfano_2009_Test_Cases';
                conj = GetAlfanoTestCases(params);
                r1(i,:) = conj.X1(1:3)';
                v1(i,:) = conj.X1(4:6)';
                C1(i,:) = [conj.C1(1,1) conj.C1(2,1) conj.C1(3,1) conj.C1(1,2) conj.C1(2,2) conj.C1(3,2) conj.C1(1,3) conj.C1(2,3) conj.C1(3,3)];
                r2(i,:) = conj.X2(1:3)';
                v2(i,:) = conj.X2(4:6)';
                C2(i,:) = [conj.C2(1,1) conj.C2(2,1) conj.C2(3,1) conj.C2(1,2) conj.C2(2,2) conj.C2(3,2) conj.C2(1,3) conj.C2(2,3) conj.C2(3,3)];
                if mod(i,3) == 1
                    % Special case 1: Zero miss distance
                    r1(i,:) = r2(i,:);
                    expSolution(i) = PcCircle(r1(i,:),v1(i,:),conj.C1,r2(i,:),v2(i,:),conj.C2,HBR(i));
                elseif mod(i,3) == 2
                    % Special case 2: Zero relative velocity
                    v1(i,:) = v2(i,:);
                    expSolution(i) = PcCircle(r1(i,:),v1(i,:),conj.C1,r2(i,:),v2(i,:),conj.C2,HBR(i));
                end
            end

            actSolution = PcCircle(r1,v1,C1,r2,v2,C2,HBR);
            for i = 1:numRows
                testCase.verifyEqual(actSolution(i),expSolution(i),'RelTol',accuracy(i));
            end
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
            
            actSolution = PcCircle(r1_ECI,v1_ECI,C1_ECI,r2_ECI,v2_ECI,C2_ECI,AlfanoHBR);
                
            testCase.verifyEqual(actSolution,AlfanoExpSolution,'RelTol',AlfanoAccuracy);
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
% E. White       | 08-09-2023 | Initial development
% L. Baars       | 12-06-2023 | Renamed to PcCircle_UnitTest.m and added
%                               some extra testing for zero miss distance
%                               and zero relative velocity testing.
% L. Baars       | 09-22-2025 | Added vectorized tests for Alfano test
%                               cases. Also, added Alfano's 12th test case
%                               (zero relative velocity) into the test set.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

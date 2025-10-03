classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Pc3D_Hall_Utils') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/OrbitTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/LoggingAndStringReporting')}) ...
        Pc3D_Hall_UnitTest < matlab.unittest.TestCase
% Pc3D_Hall_UnitTest - Unit test for Pc3D_Hall
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
% Initial version: Mar 2020;  Latest update: Jul 2024
%
% ----------------- BEGIN CODE -----------------

    properties (TestParameter)        
        AlfanoHBR = {15, 4, 15, 15, 10, 10, 10, 4, 6, 6, 4, 4};
        AlfanoTLimit = {21600, 21600, 21600, 21600, 1419, 1419, 1419, 10135, 10800, 21600, 1420, 1420};
        AlfanoAccuracy = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};
        % AlfanoAccuracy = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.005, 0.001, 0.01, 0.01, 0.01};
        % Old accuracy values, kept for reference
        AlfanoExpSolution = { 2.1680836276E-01, ...
                              1.5567748968E-02, ...
                              1.0033642320E-01, ...
                              7.3640419383E-02, ...
                              4.4489813350E-02, ...
                              4.3331309702E-03, ...
                              1.6184123020E-04, ...
                              3.5250174961E-02, ...
                              2.7983426149E-01, ...
                              3.6406379645E-01, ...
                              2.4407488450E-03, ...
                              2.4415473580E-03 };
        AlfanoFileNum = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'};
    end

    methods (Test)        
        function testOutputStruct(testCase) % Alfano Test Case 7
            HBR             = 10;
            Accuracy        = 0.005; % Desired MC Accuracy (0.01=1%)
            expNcResult     = 0.000161462; % Alfano 1E8 MC trials
            expSolution     = 0.000161462;
            
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase07.cdm');
            
            % Calculate 3D PC            
            [actSolution,outputStruct] = Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(outputStruct.Nc,expNcResult,'RelTol',Accuracy);
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
            
        end

        function testBadCovariance(testCase) % Original Test Omitron Test Case 01 (Written by Dragan Plakalovic) 
            r1      = [378.39559 4305.721887 5752.767554]*1000;
            v1      = [2.360800244 5.580331936 -4.322349039]*1000;
            r2      = [374.5180598 4307.560983 5751.130418]*1000;
            v2      = [-5.388125081 -3.946827739 3.322820358]*1000;
            cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
                     81.6751751052616  158.453402956163  -128.616921644857;
                     -67.8687662707124 -128.616921644858 105.490542562701]*1E6;
            cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
                     1.69905293875632  1.24957388457206  -1.04174164279599;
                     -1.4170164577661  -1.04174164279599 0.869260558223714]*1E6;
            HBR     = 20;
            
            % Calculate 3D PC            
            
            testCase.verifyError(@()Pc3D_Hall(r1,v1,cov1,r2,v2,cov2,HBR),...
                    'Pc3D_Hall:badCovariance')
        end

        function testNonPositiveDefinite(testCase) % Non-Pos Definite test Case
            HBR             = 20;
            Accuracy        = 1E-10; 
            expSolution     = 0;
            expNcResult     = 0;
            
           [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/OmitronTestCase_Test07_NonPDCovariance.cdm'); 
                
            % Check for non- Positive Definite Error Throw
            testCase.verifyWarning(@()PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,[],3),...
                    'RemediateCovariance2x2:remediatedNPD');
                
            % Calculate 3D PC            
            [Pc,outputStruct] = Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
            testCase.verifyEqual(Pc,expSolution,'AbsTol',Accuracy);
            testCase.verifyEqual(outputStruct.Nc,expNcResult,'AbsTol',Accuracy);
        end
        
        function test2D3DDisagreementCase(testCase) % Operational Case where 3D Pc is ~15 orders of magnitude larger than 2D estimation
            HBR             = 20;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            exp2DSolution   = 2.2660816e-20;
            exp3DSolution   = 1.9271201e-05;

            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/OmitronTestCase_Test08_3DNc.cdm'); 
            
            % Calculate 2D PC            
            [actSolution,~] = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            testCase.verifyEqual(actSolution,exp2DSolution,'RelTol',Accuracy);
            
            % Calculate 3D PC            
            [actSolution,~] = Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
            testCase.verifyEqual(actSolution,exp3DSolution,'RelTol',Accuracy);
        end
        
        % Returns a value very similar to PcElrod (see testSimilarTo3D in
        % PcElrod_UnitTest for more)
        function testSimilarTo2D(testCase)
            expSolution = 1.0289973e-02;
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
            
            [actSolution] = Pc3D_Hall(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        
        % Tests covariance matrices that are not 6x6
        function testInvalidCovarianceDimensions (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase05.cdm');
            
            % Invalid primary covariance
            testCase.verifyError(@()Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,1,r2_J2K*1000,v2_J2K*1000,C2_J2K,10),...
                    'Pc3D_Hall:badCovariance');
            % Invalid secondary covariance
            testCase.verifyError(@()Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,2,10),...
                    'Pc3D_Hall:badCovariance');
        end
        
        % Tests negative and incorrect number of HBRs
        function testInvalidHBRs (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase05.cdm');
            
            % Invalid primary HBR
            testCase.verifyError(@()Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,[-10 10]),...
                    'Pc3D_Hall:InvalidHBR');
            % Invalid secondary HBR
            testCase.verifyError(@()Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,[10 -10]),...
                    'Pc3D_Hall:InvalidHBR');
            % Invalid combined HBR
            testCase.verifyError(@()Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,-10),...
                    'Pc3D_Hall:InvalidHBR');
            % Incorrect number of HBRs
            testCase.verifyError(@()Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,[5 10 15]),...
                    'Pc3D_Hall:InvalidHBR');
        end
        
        % Tests default initialization of params
        function testDefaultParams (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase05.cdm');
            
            params = struct();
            params.gamma = 1e-16;
            params.Texpand = [];
            params.remediate_NPD_TCA_eq_covariances = false;
            params.apply_covXcorr_corrections = true;
            params.Torblimit = true;
            params.Neph = 101;
            params.POPmaxiter = 100;
            params.use_Lebedev = true;
            params.deg_Lebedev = 5810;
            params.vec_Lebedev = [];
            params.wgt_Lebedev = [];
            params.slow_method = false;
            params.AbsTol = 0;
            params.RelTol = 1e-9;
            params.MaxFunEvals = 10000;
            params.Fclip = 1e-4;
            params.Pc_tiny = 1e-300;
            params.GM = 3.986004418e14;
            params.verbose = false;
            
            [~, actOut] = Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,10,params);
            [~, expOut] = Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,10);
            
            testCase.verifyEqual(actOut.params, expOut.params);
        end
        
        % Tests invalid Texpand parameter
        function testInvalidTExpand (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase05.cdm');
            params = struct('Texpand', [1 1 1]);
            testCase.verifyError(@() Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,10,params), 'Pc3D_Hall:InvalidTExpand');
        end
        
        % Tests invalid time bounds and orbits
        function testInvalidTimeBoundsAndOrbits (testCase)
            [r1_J2K, v1_J2K, C1_J2K, r2_J2K, ~, C2_J2K] = Pc_TestSetup('../../../DataFiles/SampleCDMs/AlfanoTestCase05.cdm');
            nanV = [NaN, NaN, NaN];
            imagV = [1i, 1i, 1i];
            
            % Zero relative velocity
            testCase.verifyWarning(@() Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v1_J2K*1000,C2_J2K,10), 'Pc3D_Hall:InvalidTimeBounds');
            
            % NaN-valued velocity
            testCase.verifyWarning(@() Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,nanV,C2_J2K,10), 'EquinoctialMatrices:CalculationFailure');
            
            % Non-real-valued velocity
            testCase.verifyWarning(@() Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,imagV,C2_J2K,10), 'EquinoctialMatrices:CalculationFailure');
        end
        
        % Tests covariance cross-correlation corrections
        function testCovXCorr (testCase)
            r1 = [-2622716.15829145    1698990.97254698    7053401.88278209];
            v1 = [-4128.50471893557   -5878.79564655151    -118.60127274469];
            cov1 = [ 30945.8960599633     44172.992380352      1094.127029015    -16.9384961044964     11.1386598710759     45.9406263498005
                    44172.992380352    63081.0350093601    1600.54857635509    -24.1640670703954     15.9361771338151     65.5984249372085
                     1094.127029015    1600.54857635509    110.587616642047   -0.568332213430976    0.456096601496478     1.65805084742508
                  -16.9384961044964   -24.1640670703954  -0.568332213430976  0.00928784766563245 -0.00607332603626678  -0.0251311636623764
                   11.1386598710759    15.9361771338151   0.456096601496478 -0.00607332603626678  0.00406599298830971   0.0165654694968024
                   45.9406263498005    65.5984249372085    1.65805084742508  -0.0251311636623764   0.0165654694968024   0.0682232639995689];
            r2 = [-2631783.12714667    1685565.41148418    7053328.56776558];
            v2 = [-4112.58245687929   -5889.50355090998   -126.881531658755];
            cov2 = [ 37176528.2096518    51779295.1953059    1237287.58668701   -20025.8397475978      12497.8621510473      54618.847730921
                    51779295.195306    72120313.4426337    1723250.32582261   -27893.1756367764      17408.3723881467     76074.0995332777
                   1237287.58668701    1723250.32582261    41984.4552410924   -666.613016373449      416.862049452506      1817.7275694543
                  -20025.8397475978   -27893.1756367764   -666.613016373449    10.7880803834256     -6.73318789560721     -29.422096492432
                   12497.8621510473    17408.3723881467    416.862049452506   -6.73318789560721      4.20345007734003     18.3622071009157
                    54618.847730921    76074.0995332777     1817.7275694543    -29.422096492432      18.3622071009157     80.2453907683188];
            HBR = 14.8;
            params = struct();
            params.covXcorr.sigp = 0.116397825;
            params.covXcorr.sigs = 0.116397825;
            params.covXcorr.Gvecp = [-107.832111185126 -153.869315382605 -3.54864510814339 0.0591561311396439 -0.0385790219607754 -0.160028201582104];
            params.covXcorr.Gvecs = [-82.6219909216031 -118.567227795798 -2.89446191542217 0.0456592322551692 -0.0294390852945439 -0.123085900964496];
            
            % Covariance cross-correlation corrections applied
            [~, outTrue] = Pc3D_Hall(r1,v1,cov1,r2,v2,cov2,HBR,params);
            
            % Covariance cross-correlation corrections not applied
            [~, outFalse] = Pc3D_Hall(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyTrue(outTrue.covXcorr_corrections_applied);
            testCase.verifyFalse(outFalse.covXcorr_corrections_applied);
        end
        
    end
    
    methods (Test, ParameterCombination = 'sequential')
        
        function testAlfano(testCase, AlfanoHBR, AlfanoTLimit, AlfanoAccuracy, AlfanoExpSolution, AlfanoFileNum)
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
            
            % Bounds for 3D-Nc calculation
            TLimit = AlfanoTLimit;
            
            Nc3DParams.Neph = 1000;
            Nc3DParams.Tmin_initial = -TLimit;
            Nc3DParams.Tmax_initial =  TLimit;
            Nc3DParams.Tmin_limit   = Nc3DParams.Tmin_initial;
            Nc3DParams.Tmax_limit   = Nc3DParams.Tmax_initial;
            
            actSolution = Pc3D_Hall(r1_ECI,v1_ECI,C1_ECI,r2_ECI,v2_ECI,C2_ECI,AlfanoHBR,Nc3DParams);
            
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
% T. Lechtenberg | 03-16-2020 | Initial Development
% T. Lechtenberg | 06-09-2020 | Added unit test representing operational
%                               event where 2D Pc was 15 orders of
%                               magnitude less than 3D Pc
% T. Lechtenberg | 08-03-2020 | Changed errors caught in some situations
%                               and changed NPD expRes from NaN to 0
% L. Baars       | 10-03-2022 | Changed Alfano test case file reading
%                               changed path fixtures
% L. Baars       | 02-22-2023 | Changed path fixtures, fixed unit test
%                               cases after Pc code reorganization
% E. White       | 08-08-2023 | Added compliant documentation, greatly
%                               improved readability and decreased length
%                               of unit tests, changed unit tests to use
%                               correct Alfano data
% L. Baars       | 12-27-2023 | Updated unit test to reflect changes of
%                               errors to warnings within Pc3D_Hall.m
% L. Baars       | 07-23-2024 | Fixed warning code for one of the checks.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

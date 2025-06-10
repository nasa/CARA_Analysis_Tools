classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../../AugmentedMath') ...
        matlab.unittest.fixtures.PathFixture('../../CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../TimeTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../CovarianceTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../PosVelTransformations') ...
        matlab.unittest.fixtures.PathFixture('../../../ProbabilityOfCollision') ...
        matlab.unittest.fixtures.PathFixture('../../../ProbabilityOfCollision/Utils')}) ...
        Pc_CovXCorr_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the Determination of TCA rectification from if improperly characterized in state estimates
    methods (Test)
        function test01(testCase) % Non-Zero De-Correlated Pc
            Accuracy        = 5E-3;
            expXCorrPc      = 8.16E-06;
            HBR             = 6;
            % ID of Test Event
            FileName = 'InputFiles/DensityDecorrelationTestCaseCDM01.txt';
            
            % Parse CDM Input File
            [cdmhead, cdmobj] = read_cdm(FileName);

            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end
            
            % Set up Primary Object State and Covariance
            x1         = [cdmobj(1).X cdmobj(1).Y cdmobj(1).Z cdmobj(1).X_DOT cdmobj(1).Y_DOT cdmobj(1).Z_DOT]*1000;
            % Convert Primary Object State to Inertial Frame if required
            if strcmpi(cdmobj(1).REF_FRAME,'ITRF')
                [r1_J2K,v1_J2K] = PosVelConvert(x1(1:3)/1000,x1(4:6)/1000,strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
            elseif strcmpi(cdmobj(1).REF_FRAME,'EME2000')
                r1_J2K = x1(1:3)/1000;
                v1_J2K = x1(4:6)/1000;
            else
                error('Incorrect or Non-Existent Reference Frame Specified for Primary Object State')
            end

            % Add Matrices to CDM Object
            CDM.OBJECT1_R_J2K = r1_J2K;
            CDM.OBJECT1_V_J2K = v1_J2K;

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

            % Add Matrices to CDM Object
            CDM.OBJECT1_C_RTN = C1_RTN;
            CDM.OBJECT1_C_J2K = C1_J2K;

            % Set up Secondary Object State and Covariance
            x2         = [cdmobj(2).X cdmobj(2).Y cdmobj(2).Z cdmobj(2).X_DOT cdmobj(2).Y_DOT cdmobj(2).Z_DOT]*1000;
            % Convert Secondary Object State to Inertial Frame
            if strcmpi(cdmobj(1).REF_FRAME,'ITRF')
                [r2_J2K,v2_J2K] = PosVelConvert(x2(1:3)/1000,x2(4:6)/1000,strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
            elseif strcmpi(cdmobj(1).REF_FRAME,'EME2000')
                r2_J2K = x2(1:3)/1000;
                v2_J2K = x2(4:6)/1000;
            else
                error('Incorrect or Non-Existent Reference Frame Specified for Secondary Object State')
            end

            % Add Matrices to CDM Object
            CDM.OBJECT2_R_J2K = r2_J2K;
            CDM.OBJECT2_V_J2K = v2_J2K;

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K);

            % Add Matrices to CDM Object
            CDM.OBJECT2_C_RTN = C2_RTN;
            CDM.OBJECT2_C_J2K = C2_J2K;
            
            % Calculate 0 Correlation Pc           
            [Pc_Elrod,~,~,~] = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
            
            % Calculate Decorrelated Pc
            [PcXC,CovXC,~]  = Pc_CovXCorr(cdmhead,cdmobj,HBR,1);
                
            % Verify Expected Solution
            testCase.verifyEqual(PcXC,expXCorrPc,'RelTol',Accuracy);
            
            % Verify Pc Estimate Using Output Covariance Matches Output
            % Covariance
            [Pc_ElrodXC,~,~,~] = PcElrod(r1_J2K*1000,v1_J2K*1000,CovXC,r2_J2K*1000,v2_J2K*1000,zeros(3,3),HBR);
            testCase.verifyEqual(PcXC,Pc_ElrodXC,'RelTol',Accuracy);
            
            % Verify Pc Estimate with 0(%) density consider parameter
            % matches output without
            cdmobj(1).DCP_DENSITY_UNCERTAINTY = 0;
            cdmobj(2).DCP_DENSITY_UNCERTAINTY = 0;
            [PcXC,CovXC,~]  = Pc_CovXCorr(cdmhead,cdmobj,HBR,1);
            testCase.verifyEqual(PcXC,Pc_Elrod,'RelTol',Accuracy);
            
        end
        
        % Commented out as logic to determine this is contained within
        % "PC_Omnibus.m"
%         function test02(testCase) % CDM with only one object with decorrelation data
%             Accuracy        = 5E-3;
%             expXCorrPc      = 4.60E-05;
%             HBR             = 6;
%             % ID of Test Event
%             FileName = 'InputFiles/DensityDecorrelationTestCaseCDM02.txt';
%             
%             % Parse CDM Input File
%             [cdmhead, cdmobj] = read_cdm(FileName);
% 
%             % Process CDM head
%             % Get available FieldNames from CDMhead object
%             names                       = fieldnames(cdmhead);
%             for j=1:length(names)
%                 CDM.(names{j})   = cdmhead.(names{j});
%             end
% 
%             % Process CDM Object Entries
%             % Get available FieldNames from CDMobject
%             names                       = fieldnames(cdmobj);
%             for j=1:length(names)
%                 CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
%             end
%             for j=1:length(names)
%                 CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
%             end
%             
%             % Set up Primary Object State and Covariance
%             x1         = [cdmobj(1).X cdmobj(1).Y cdmobj(1).Z cdmobj(1).X_DOT cdmobj(1).Y_DOT cdmobj(1).Z_DOT]*1000;
%             % Convert Primary Object State to Inertial Frame if required
%             if strcmpi(cdmobj(1).REF_FRAME,'ITRF')
%                 [r1_J2K,v1_J2K] = PosVelConvert(x1(1:3)/1000,x1(4:6)/1000,strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
%             elseif strcmpi(cdmobj(1).REF_FRAME,'EME2000')
%                 r1_J2K = x1(1:3)/1000;
%                 v1_J2K = x1(4:6)/1000;
%             else
%                 error('Incorrect or Non-Existent Reference Frame Specified for Primary Object State')
%             end
% 
%             % Add Matrices to CDM Object
%             CDM.OBJECT1_R_J2K = r1_J2K;
%             CDM.OBJECT1_V_J2K = v1_J2K;
% 
%             % Get Primary Object RTN Covariance
%             C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
%                        cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
%                        cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
%                        cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
%                        cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
%                        cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
%             % Transform Primary Object Covariance to J2K
%             [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);
% 
%             % Add Matrices to CDM Object
%             CDM.OBJECT1_C_RTN = C1_RTN;
%             CDM.OBJECT1_C_J2K = C1_J2K;
% 
%             % Set up Secondary Object State and Covariance
%             x2         = [cdmobj(2).X cdmobj(2).Y cdmobj(2).Z cdmobj(2).X_DOT cdmobj(2).Y_DOT cdmobj(2).Z_DOT]*1000;
%             % Convert Secondary Object State to Inertial Frame
%             if strcmpi(cdmobj(1).REF_FRAME,'ITRF')
%                 [r2_J2K,v2_J2K] = PosVelConvert(x2(1:3)/1000,x2(4:6)/1000,strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
%             elseif strcmpi(cdmobj(1).REF_FRAME,'EME2000')
%                 r2_J2K = x2(1:3)/1000;
%                 v2_J2K = x2(4:6)/1000;
%             else
%                 error('Incorrect or Non-Existent Reference Frame Specified for Secondary Object State')
%             end
% 
%             % Add Matrices to CDM Object
%             CDM.OBJECT2_R_J2K = r2_J2K;
%             CDM.OBJECT2_V_J2K = v2_J2K;
% 
%             C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
%                        cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
%                        cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
%                        cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
%                        cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
%                        cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
%             % Transform Secondary Object Covariance to J2K
%             [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K);
% 
%             % Add Matrices to CDM Object
%             CDM.OBJECT2_C_RTN = C2_RTN;
%             CDM.OBJECT2_C_J2K = C2_J2K;
%             
%             % Calculate 0 Correlation Pc           
%             [Pc_Elrod,~,~,~] = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
%             
%             % Calculate Decorrelated Pc
%             [PcXC,CovXC,~]  = Pc_CovXCorr(cdmhead,cdmobj,HBR,1);
%                 
%             % Verify Expected Solution
%             testCase.verifyEqual(PcXC,expXCorrPc,'RelTol',Accuracy);
%             
%             % Verify Pc Estimate Using Output Covariance Matches Output
%             % Covariance
%             [Pc_ElrodXC,~,~,~] = PcElrod(r1_J2K*1000,v1_J2K*1000,CovXC,r2_J2K*1000,v2_J2K*1000,zeros(3,3),HBR);
%             testCase.verifyEqual(PcXC,Pc_ElrodXC,'RelTol',Accuracy);
%             
%             % Verify Pc Estimate with 0(%) density consider parameter
%             % matches output without
%             cdmobj(1).DCP_DENSITY_UNCERTAINTY = 0;
%             cdmobj(2).DCP_DENSITY_UNCERTAINTY = 0;
%             [PcXC,CovXC,~]  = Pc_CovXCorr(cdmhead,cdmobj,HBR,1);
%             testCase.verifyEqual(PcXC,Pc_Elrod,'RelTol',Accuracy);
%             
%         end

        function test03(testCase) % Similar to Event in Test 02 but both objects have Decorrelation Parameters
            Accuracy        = 5E-3;
            expXCorrPc      = 3.20E-05;
            HBR             = 6;
            % ID of Test Event
            FileName = 'InputFiles/DensityDecorrelationTestCaseCDM03.txt';
            
            % Parse CDM Input File
            [cdmhead, cdmobj] = read_cdm(FileName);

            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end
            
            % Set up Primary Object State and Covariance
            x1         = [cdmobj(1).X cdmobj(1).Y cdmobj(1).Z cdmobj(1).X_DOT cdmobj(1).Y_DOT cdmobj(1).Z_DOT]*1000;
            % Convert Primary Object State to Inertial Frame if required
            if strcmpi(cdmobj(1).REF_FRAME,'ITRF')
                [r1_J2K,v1_J2K] = PosVelConvert(x1(1:3)/1000,x1(4:6)/1000,strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
            elseif strcmpi(cdmobj(1).REF_FRAME,'EME2000')
                r1_J2K = x1(1:3)/1000;
                v1_J2K = x1(4:6)/1000;
            else
                error('Incorrect or Non-Existent Reference Frame Specified for Primary Object State')
            end

            % Add Matrices to CDM Object
            CDM.OBJECT1_R_J2K = r1_J2K;
            CDM.OBJECT1_V_J2K = v1_J2K;

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

            % Add Matrices to CDM Object
            CDM.OBJECT1_C_RTN = C1_RTN;
            CDM.OBJECT1_C_J2K = C1_J2K;

            % Set up Secondary Object State and Covariance
            x2         = [cdmobj(2).X cdmobj(2).Y cdmobj(2).Z cdmobj(2).X_DOT cdmobj(2).Y_DOT cdmobj(2).Z_DOT]*1000;
            % Convert Secondary Object State to Inertial Frame
            if strcmpi(cdmobj(1).REF_FRAME,'ITRF')
                [r2_J2K,v2_J2K] = PosVelConvert(x2(1:3)/1000,x2(4:6)/1000,strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
            elseif strcmpi(cdmobj(1).REF_FRAME,'EME2000')
                r2_J2K = x2(1:3)/1000;
                v2_J2K = x2(4:6)/1000;
            else
                error('Incorrect or Non-Existent Reference Frame Specified for Secondary Object State')
            end

            % Add Matrices to CDM Object
            CDM.OBJECT2_R_J2K = r2_J2K;
            CDM.OBJECT2_V_J2K = v2_J2K;

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K);

            % Add Matrices to CDM Object
            CDM.OBJECT2_C_RTN = C2_RTN;
            CDM.OBJECT2_C_J2K = C2_J2K;
            
            % Calculate 0 Correlation Pc           
            [Pc_Elrod,~,~,~] = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
            
            % Calculate Decorrelated Pc
            [PcXC,CovXC,~]  = Pc_CovXCorr(cdmhead,cdmobj,HBR,1);
                
            % Verify Expected Solution
            testCase.verifyEqual(PcXC,expXCorrPc,'RelTol',Accuracy);
            
            % Verify Pc Estimate Using Output Covariance Matches Output
            % Covariance
            [Pc_ElrodXC,~,~,~] = PcElrod(r1_J2K*1000,v1_J2K*1000,CovXC,r2_J2K*1000,v2_J2K*1000,zeros(3,3),HBR);
            testCase.verifyEqual(PcXC,Pc_ElrodXC,'RelTol',Accuracy);
            
            % Verify Pc Estimate with 0(%) density consider parameter
            % matches output without
            cdmobj(1).DCP_DENSITY_UNCERTAINTY = 0;
            cdmobj(2).DCP_DENSITY_UNCERTAINTY = 0;
            [PcXC,CovXC,~]  = Pc_CovXCorr(cdmhead,cdmobj,HBR,1);
            testCase.verifyEqual(PcXC,Pc_Elrod,'RelTol',Accuracy);
            
        end
    end 
end
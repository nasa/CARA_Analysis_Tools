classdef PcElrod_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the 2D probability of Collision Calculation at TCA

    methods (Test)
        function test01(testCase) % Alfano Test Case 1
            HBR             = 15;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.146749549;
            % ID of Test Event
            FileName = 'AlfanoTestCase01.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test02(testCase) % Alfano Test Case 2
            HBR             = 4;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.006222267;
            % ID of Test Event
            FileName = 'AlfanoTestCase02.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test03(testCase) % Alfano Test Case 3
            HBR             = 15;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.100351176;
            % ID of Test Event
            FileName = 'AlfanoTestCase03.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test04(testCase) % Alfano Test Case 4
            HBR             = 15;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.049323406;
            % ID of Test Event
            FileName = 'AlfanoTestCase04.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,64);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test05(testCase) % Alfano Test Case 5
            HBR             = 10;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.044487386;
            % ID of Test Event
            FileName = 'AlfanoTestCase05.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,32);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test06(testCase) % Alfano Test Case 6
            HBR             = 10;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.004335455;
            % ID of Test Event
            FileName = 'AlfanoTestCase06.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test07(testCase) % Alfano Test Case 7
            HBR             = 10;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.000158147;
            % ID of Test Event
            FileName = 'AlfanoTestCase07.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test08(testCase) % Alfano Test Case 8
            HBR             = 4;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.036948008;
            % ID of Test Event
            FileName = 'AlfanoTestCase08.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test09(testCase) % Alfano Test Case 9
            HBR             = 6;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.290146291;
            % ID of Test Event
            FileName = 'AlfanoTestCase09.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,48);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test10(testCase) % Alfano Test Case 10
            HBR             = 6;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.290146291;
            % ID of Test Event
            FileName = 'AlfanoTestCase10.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,48);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test11(testCase) % Alfano Test Case 11
            HBR             = 4;
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expSolution     = 0.002672026;
            % ID of Test Event
            FileName = 'AlfanoTestCase11.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
            
            % Calculate 2D PC            
            actSolution = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test12(testCase) % Original Test Omitron Test Case 01 (Written by Dragan Plakalovic)
            % Circular Area
            expSolution      = 2.70601573490125e-05;
            Accuracy        = 0.001; 
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
            Tol     = 1e-09;
            HBRType = 'circle';
            [actSolution]    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
%  %%%%%%Following cases not used because PcGaussQuad error function
%  integration doesn't accomodate square areas.
%       function test13(testCase) % Original Test Omitron Test Case 02 (Written by Dragan Plakalovic)
%             % Square Area
%             expSolution      = 3.44534649703584e-05;
%             Accuracy        = 0.001; 
%             r1      = [378.39559 4305.721887 5752.767554];
%             v1      = [2.360800244 5.580331936 -4.322349039];
%             r2      = [374.5180598 4307.560983 5751.130418];
%             v2      = [-5.388125081 -3.946827739 3.322820358];
%             cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
%                      81.6751751052616  158.453402956163  -128.616921644857;
%                      -67.8687662707124 -128.616921644858 105.490542562701];
%             cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
%                      1.69905293875632  1.24957388457206  -1.04174164279599;
%                      -1.4170164577661  -1.04174164279599 0.869260558223714];
%             HBR     = 0.020;
%             Tol     = 1e-09;
%             HBRType = 'square';
%             [actSolution]    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR,Tol,HBRType);
%             
%             testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
%         end
%         function test14(testCase) % Original Test Omitron Test Case 03 (Written by Dragan Plakalovic)
%             % Square Eqiuvalent Area
%             expSolution      = 2.70601573490111e-05;
%             Accuracy        = 0.001; 
%             r1      = [378.39559 4305.721887 5752.767554];
%             v1      = [2.360800244 5.580331936 -4.322349039];
%             r2      = [374.5180598 4307.560983 5751.130418];
%             v2      = [-5.388125081 -3.946827739 3.322820358];
%             cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
%                      81.6751751052616  158.453402956163  -128.616921644857;
%                      -67.8687662707124 -128.616921644858 105.490542562701];
%             cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
%                      1.69905293875632  1.24957388457206  -1.04174164279599;
%                      -1.4170164577661  -1.04174164279599 0.869260558223714];
%             HBR     = 0.020;
%             Tol     = 1e-09;
%             HBRType = 'squareEquArea';
%             [actSolution]    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR,Tol,HBRType);
%             
%             testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
%         end
        function test13(testCase) % Non-Pos Definite test Case
            HBR             = 20;
            Accuracy        = 1E-10; 
            expSolution     = 0;
            % ID of Test Event
            FileName = 'OmitronTestCase_Test07_NonPDCovariance.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            [cdmhead, cdmobj] = read_cdm(fullfile(PathName, FileName));
            
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

            % Get Primary Object RTN Covariance
            C1_RTN  = [cdmobj(1).CR_R     cdmobj(1).CT_R    cdmobj(1).CN_R    cdmobj(1).CRDOT_R    cdmobj(1).CTDOT_R    cdmobj(1).CNDOT_R
                       cdmobj(1).CT_R     cdmobj(1).CT_T    cdmobj(1).CN_T    cdmobj(1).CRDOT_T    cdmobj(1).CTDOT_T    cdmobj(1).CNDOT_T
                       cdmobj(1).CN_R     cdmobj(1).CN_T    cdmobj(1).CN_N    cdmobj(1).CRDOT_N    cdmobj(1).CTDOT_N    cdmobj(1).CNDOT_N
                       cdmobj(1).CRDOT_R  cdmobj(1).CRDOT_T cdmobj(1).CRDOT_N cdmobj(1).CRDOT_RDOT cdmobj(1).CTDOT_RDOT cdmobj(1).CNDOT_RDOT
                       cdmobj(1).CTDOT_R  cdmobj(1).CTDOT_T cdmobj(1).CTDOT_N cdmobj(1).CTDOT_RDOT cdmobj(1).CTDOT_TDOT cdmobj(1).CNDOT_TDOT
                       cdmobj(1).CNDOT_R  cdmobj(1).CNDOT_T cdmobj(1).CNDOT_N cdmobj(1).CNDOT_RDOT cdmobj(1).CNDOT_TDOT cdmobj(1).CNDOT_NDOT];
            % Transform Primary Object Covariance to J2K
            [C1_J2K] = RIC2ECI(C1_RTN,r1_J2K,v1_J2K);

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

            C2_RTN  = [cdmobj(2).CR_R     cdmobj(2).CT_R    cdmobj(2).CN_R    cdmobj(2).CRDOT_R    cdmobj(2).CTDOT_R    cdmobj(2).CNDOT_R
                       cdmobj(2).CT_R     cdmobj(2).CT_T    cdmobj(2).CN_T    cdmobj(2).CRDOT_T    cdmobj(2).CTDOT_T    cdmobj(2).CNDOT_T
                       cdmobj(2).CN_R     cdmobj(2).CN_T    cdmobj(2).CN_N    cdmobj(2).CRDOT_N    cdmobj(2).CTDOT_N    cdmobj(2).CNDOT_N
                       cdmobj(2).CRDOT_R  cdmobj(2).CRDOT_T cdmobj(2).CRDOT_N cdmobj(2).CRDOT_RDOT cdmobj(2).CTDOT_RDOT cdmobj(2).CNDOT_RDOT
                       cdmobj(2).CTDOT_R  cdmobj(2).CTDOT_T cdmobj(2).CTDOT_N cdmobj(2).CTDOT_RDOT cdmobj(2).CTDOT_TDOT cdmobj(2).CNDOT_TDOT
                       cdmobj(2).CNDOT_R  cdmobj(2).CNDOT_T cdmobj(2).CNDOT_N cdmobj(2).CNDOT_RDOT cdmobj(2).CNDOT_TDOT cdmobj(2).CNDOT_NDOT];
            % Transform Secondary Object Covariance to J2K
            [C2_J2K] = RIC2ECI(C2_RTN,r2_J2K,v2_J2K); 
                
            % Check for non- Positive Definite Error Throw
            warning on
            testCase.verifyWarning(@()PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR),...
                    'PcElrod:accNPD')
            warning off
            [Pc,~,IsPosDef,IsRemediated] = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
            testCase.verifyEqual(Pc,expSolution,'AbsTol',Accuracy);
            testCase.verifyTrue(IsPosDef);
            testCase.verifyTrue(IsRemediated);
        end
    end 
end
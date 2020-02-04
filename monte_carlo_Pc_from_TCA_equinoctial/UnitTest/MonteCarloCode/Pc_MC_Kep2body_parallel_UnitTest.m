classdef Pc_MC_Kep2body_parallel_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the rectilinear Monte Carlo Sampling at TCA 
    
    methods (Test)
        function test01(testCase) % High 2D Probability of Collision
            HBR             = 20;
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 4.20E-01;
            % ID of Test Event
            FileName = 'OmitronTestCase_Test01_HighPc.cdm';
            
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test02(testCase) % High Radial Position Uncertainty
            HBR             = 100; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 3.0216E-03;
            % ID of Test Event
            FileName = 'OmitronTestCase_Test02_MaxRadialSigma.cdm';
            
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test03(testCase) % High Intrack Position Uncertainty
            HBR             = 100; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 2.6516E-03;
            % ID of Test Event
            FileName = 'OmitronTestCase_Test03_MaxIntrackSigma.cdm';
            
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test04(testCase) % High Crosstrack Position Uncertainty
            HBR             = 100; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 2.4989E-03;
            % ID of Test Event
            FileName = 'OmitronTestCase_Test04_MaxCrossTrackSigma.cdm';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            CDMfile  = fullfile(PathName,FileName);
            
            if exist(CDMfile,'file') % CDM file removed from SDKs for public distribution, returns a test pass if CDM file not present
                [cdmhead, cdmobj] = read_cdm(CDMfile);

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
                Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');

                % Estimate NUmer of Samples
                Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
                Nsample_kep = max(1e2,Nsample_kep);
                Nsample_kep = min(1e10,Nsample_kep);

                % Determine Batch Size min 1000, max 5000
                p = gcp;
                Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);

                % Get Primary Objecty Orbital Period
                a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
                OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

                % Get Bounds
                [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);

                % Check if event has low relative velocity
                if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                    warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                    warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                    tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                    tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                    fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
                end
                tfc_bins.N   = 100;
                tfc_bins.x1  = tau0;
                tfc_bins.x2  = tau1;
                tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
                tfc_bins.del = tfc_bins.wid/tfc_bins.N;
                tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
                tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
                tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

                % Run MC Code
                X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
                [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
                if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                    fr      = -1;
                else
                    fr      = 1;
                end
                [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
                E1TCA   = [n af ag chi psi lM fr];
                Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
                Jctoe1 = inv(Jctoe1);
                PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

                X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
                [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
                if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                    fr      = -1;
                else
                    fr      = 1;
                end
                [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
                E2TCA   = [n af ag chi psi lM fr];
                Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
                Jctoe2 = inv(Jctoe2);
                PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA

                % Equinoctial Sampling
                [~, ~, ~,...
                     actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                     ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                        Nsample_kep, tfc_bins,...
                        E1TCA', [], PEq1TCA,...
                        E2TCA', [], PEq2TCA,...
                        HBR, GM, 'k2bpri',...
                        MC_Confidence, Nsample_batch);

                testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
            else
                testCase.verifyEqual(1,1,'RelTol',Accuracy);
            end
        end
        function test05(testCase) % Low Miss
            HBR             = 20; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 1.70E-03;
            % ID of Test Event
            FileName = 'OmitronTestCase_Test05_MinMiss.cdm';
            
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test06(testCase) % Low Relative Velocity
            HBR             = 20; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 1.1325E-01;
            % ID of Test Event
            FileName = 'OmitronTestCase_Test06_MinRelVel.cdm';
            
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test07(testCase) % Alfano Test Case 1
            HBR             = 15; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 0.217467140;
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                tau1    = 21600;
                tau0    = -21600;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test08(testCase) % Alfano Test Case 3
            HBR             = 15; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 0.100846420;
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test09(testCase) % Alfano Test Case 4
            HBR             = 15; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 0.073089530;
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tau1    = 21600;
            tau0    = -21600;
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test10(testCase) % Alfano Test Case 5
            HBR             = 10; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 0.044498913;
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tau1    = 1420;
            tau0    = -1420;
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test11(testCase) % Alfano Test Case 6
            HBR             = 10; % Modified to result in reduced number of required Monte Carlo Trials (Pc>1e-3)
            Accuracy        = 0.1; % Desired MC Accuracy (0.01=1%)
            MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
            GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2] 
            LowVelThreshold = 0.05; % 5% of Orbital Period
            expSolution     = 0.004300500;
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
            Pc2D = Pc2D_Foster(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR,1e-8,'circle');
            
            % Estimate NUmer of Samples
            Nsample_kep = EstimateRequiredSamples(Pc2D,Accuracy*0.75,MC_Confidence);
            Nsample_kep = max(1e2,Nsample_kep);
            Nsample_kep = min(1e10,Nsample_kep);
            
            % Determine Batch Size min 1000, max 5000
            p = gcp;
            Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);
            
            % Get Primary Objecty Orbital Period
            a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
            OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

            % Get Bounds
            [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
            
            % Check if event has low relative velocity
            if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
                warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
                warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
                tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
                tau0    = tau1-LowVelThreshold*OrbitalPeriod;
                fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
            end
            tau1    = 1420;
            tau0    = -1420;
            tfc_bins.N   = 100;
            tfc_bins.x1  = tau0;
            tfc_bins.x2  = tau1;
            tfc_bins.wid = (tfc_bins.x2-tfc_bins.x1);
            tfc_bins.del = tfc_bins.wid/tfc_bins.N;
            tfc_bins.xhi = tfc_bins.x1+(1:tfc_bins.N)*tfc_bins.del;
            tfc_bins.xlo = tfc_bins.xhi-tfc_bins.del;
            tfc_bins.xmd = tfc_bins.xhi-tfc_bins.del/2;

            % Run MC Code
            X1TCA   = [r1_J2K v1_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr);
            E1TCA   = [n af ag chi psi lM fr];
            Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr);
            Jctoe1 = inv(Jctoe1);
            PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

            X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
            [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad');
            if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
                fr      = -1;
            else
                fr      = 1;
            end
            [~,n,af,ag,chi,psi,lM,F] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr);
            E2TCA   = [n af ag chi psi lM fr];
            Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr);
            Jctoe2 = inv(Jctoe2);
            PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA
            
            % Equinoctial Sampling
            [~, ~, ~,...
                 actSolution, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
                 ~  , ~  , ~] = Pc_MC_Kep2body_parallel(...
                    Nsample_kep, tfc_bins,...
                    E1TCA', [], PEq1TCA,...
                    E2TCA', [], PEq2TCA,...
                    HBR, GM, 'k2bpri',...
                    MC_Confidence, Nsample_batch);
                
            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
    end
    
end
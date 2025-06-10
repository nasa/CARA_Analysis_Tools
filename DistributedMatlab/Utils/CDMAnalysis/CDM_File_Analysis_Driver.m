%% CDM File Analysis Driver
% 
% Purpose: This code is intended to act as a stand alone Driver to parse and process
% an input CDM txt file and report Probability of Collision estimates using
% multiple Methods

%% Set Up File Paths
p = mfilename('fullpath');
[filepath,~,~] = fileparts(p);
cd([filepath filesep '..' filesep '..']);
addpath(genpath(pwd))
ParentPath = pwd;

%% Constants
MC_Confidence   = 0.95; % Desired Monte Carlo Confidence Level
GM              = 3.986004418e14;% Earth gravitational constant mu = GM (EGM-96) [m^3/s^2]
LowVelThreshold = 0.05; % 5% of Orbital Period
DefaultCd       = 2.7; % Default Drag Coefficient Used for Estimating Masses of Primary and Secondary Objects
DefaultCdSigma  = 0.10; % Relative Cd uncertainty for use in Mass Estimation Process
MassQuantiles   = [0.50 0.68 0.75 0.80 0.90 0.95 0.99 0.999 0.9999]; % Mass Estimation Quantiles
DefaultQuantile = 0.999; % Default Secondary Mass Estimation Quantile for Use in Collision Consequence
MassEstimateFile= fullfile(ParentPath,'Main','DataFiles','MassQuantileEstimates(kg).txt');

ObjectTypeTable = [{'PAYLOAD'}          1
                   {'ROCKET BODY'}      2
                   {'DEBRIS'}           3
                   {'UNKNOWN'}          4
                   {'OTHER'}            5
                   {'INACTIVE PAYLOAD'} 6
                   {'ACTIVE PAYLOAD'}   7];
                   
ODQAFlags = {'Geopotential Model' 'Proper Drag Application' 'Proper SRP Application'...
             'Lunar/Solar Perturbation Application' 'Earth Tides Application' ...
             'Unreasonable Drag Coefficient' 'Unreasonable SRP Coefficient'...
             'LUPI Ratio' 'Residual Acceptance' 'WRMS Value' 'ASW Default Covariance'...
             'Covariance Matrix Not Positive Definite' 'Intrack Position Uncertainty'...
             'Epoch Age'};

%% Allow User to Select Which Services They Want to Perform
AvailableServices = {'Pc Omnibus Tool',                     true;
                     'Single Covariance Maximum 2D Pc',     true;
                     'Secondary Object Mass Estimation',    true;
                     'Collision Consequence',               true;
                     'CDM Secondary Object ODQA',           true;
                     'Monte Carlo Probability of Collision',false};

AvailableServices = ServiceSelectionGUI(AvailableServices);
AvailableServices = cell2struct(AvailableServices(:,2),strrep(AvailableServices(:,1),' ','')');
%% Allow User to Select CDM and Set Inputs
% Retrieve File
[FileName,PathName] = uigetfile({'*.txt;*.cdm';'*.*'},'Please Select A CDM Input File');
if FileName==0
    fprintf('\nNo file selected, program will now exit\n')
    return
end

% Parse CDM Input File
CDMfile = fullfile(PathName, FileName);
[cdmhead, cdmobj] = read_cdm(CDMfile);

% Get Hard Body Radius (m) if required
if (isfield(AvailableServices,'PcOmnibusTool') && AvailableServices.PcOmnibusTool)||...
   (isfield(AvailableServices,'SingleCovarianceMaximum2DPc') && AvailableServices.SingleCovarianceMaximum2DPc) ||...
   (isfield(AvailableServices,'MonteCarloProbabilityofCollision') && AvailableServices.MonteCarloProbabilityofCollision)
        HBR         = 20;
        answer      = inputdlg('Input Hard Body Radius (HBR) in Meters','HBR',1,{num2str(HBR)});
        if isempty(answer)
            fprintf('\nNo Hard Body Radius Input, program will now exit\n')
            return
        end
        HBR         = str2double(answer{1});
end

if isfield(AvailableServices,'CollisionConsequence') && AvailableServices.CollisionConsequence
    % Estimate and request confirmation of primary object mass using a default Cd 0 if required
    [PrimaryMass,~] = EstimateMassFromRCS(cdmobj(1).AREA_PC,DefaultCd,cdmobj(1).CD_AREA_OVER_MASS);
    answer          = inputdlg('Input Primary Object Mass - kg (Default estimated from drag characteristics)'...
                                ,'Primary Object Mass',1,{num2str(PrimaryMass,'%0.2f')});
    if isempty(answer)
        fprintf('\nNo Primary Mass Input, program will not estimate consequence of collision\n')
        PrimaryMass         = answer;
    else
        PrimaryMass         = str2double(answer{1});
        if isnan(PrimaryMass)
            fprintf('\nNo Primary Mass Input, program will not estimate consequence of collision\n')
            PrimaryMass         = [];
        end
    end
end

% Get Desired Monte Carlo Accuracy if required
if isfield(AvailableServices,'MonteCarloProbabilityofCollision') && AvailableServices.MonteCarloProbabilityofCollision
    Accuracy    = 0.1; % Desired MC Accuracy (0.01=1%)
    answer      = inputdlg('Input Desired Accuracy of Monte Carlo Simulation (Precision of Answer to percieved "truth" value)','Monte Carlo Accuracy',1,{[num2str(Accuracy*100,'%0.2f') '%']});
    if isempty(answer)
        fprintf('\nNo Desired Accuracy Input, program will now exit\n')
        return
    end
    answer      = strsplit(answer{1},'%');
    Accuracy    = str2double(answer{1})/100;
end

% Set Output folder for generated Figures if required
if (isfield(AvailableServices,'MonteCarloProbabilityofCollision') && AvailableServices.MonteCarloProbabilityofCollision) ||...
    (isfield(AvailableServices,'PcOmnibusTool') && AvailableServices.PcOmnibusTool)
        FiguresFolder = uigetdir('','Please Select A Directory to save figures to (optional)');
        if FiguresFolder==0
            fprintf('\nNo output folder selected, program will continue though plots will not be generated\n')
        end
end
% FiguresFolder   = 'D:\TemporaryFigures\';
%% Set up and convert primary and secondary states and covariances if PC Omnibus Tool Not Selected
if (~isfield(AvailableServices,'PcOmnibusTool') || ~AvailableServices.PcOmnibusTool)
    
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
else
%% Execute Pc Omnibus Tool
    [PcStructOut]   = Pc_Omnibus(CDMfile,HBR,FiguresFolder);
    
    % Populate Variables as needed
    CDM             = PcStructOut.CDM;
    r1_J2K          = CDM.OBJECT1_R_J2K;
    v1_J2K          = CDM.OBJECT1_V_J2K;
    C1_J2K          = CDM.OBJECT1_C_J2K;
    r2_J2K          = CDM.OBJECT2_R_J2K;
    v2_J2K          = CDM.OBJECT2_V_J2K;
    C2_J2K          = CDM.OBJECT2_C_J2K;
    Pc_Elrod        = PcStructOut.PC2DBase;
end
%% Calculate Maximum Possible 2D PC 
if (isfield(AvailableServices,'SingleCovarianceMaximum2DPc') && AvailableServices.SingleCovarianceMaximum2DPc)
    MaxPc = FrisbeeMaxPc(r1_J2K*1000,v1_J2K*1000,C1_J2K(1:3,1:3),r2_J2K*1000,v2_J2K*1000,C2_J2K(1:3,1:3),HBR,1e-8,'circle');

    % Report Calculated 2D Pc
    fprintf(['\nMaximum Possible 2D probability of collision (Frisbee Formulation) from input \n\tCDM assuming one object covariance is unpopulated: \t' num2str(MaxPc,'%.4E') '\n\n'])
end
%% Calculate Secondary Mass Estimates
clear QuantileArray
if isfield(AvailableServices,'SecondaryObjectMassEstimation') && AvailableServices.SecondaryObjectMassEstimation
    RCS         = CDM.OBJECT2_AREA_PC; 
    if RCS == 0
        warning('Secondary RCS is listed as "0" in CDM, Please review CDM file, Secondary object mass will be estimated using a minimal RCS value of: 0.0001 m^2');
        RCS = 0.0001;
    end
    B           = CDM.OBJECT2_CD_AREA_OVER_MASS;

    % Check if ballistic coefficient variance is populated
    if ~isfield(cdmobj,'CDRG_DRG') ||...
            isempty(CDM.OBJECT2_CDRG_DRG)|| ...%isNull(CDM.OBJECT2_CDRG_DRG)|| ...
            CDM.OBJECT2_CDRG_DRG == 100
        BVar        = 0;
    else
        BVar        = CDM.OBJECT2_CDRG_DRG;
    end
    
    Cd          = DefaultCd;
    CdVar       = (DefaultCd*DefaultCdSigma)^2; % Get Cd Variance from Relative Uncertainty Specified above
    [massVec,QuantileArray] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,MassQuantiles,[]);
    
    for i=1:length(QuantileArray)
        fprintf(['Secondary Mass Estimate at ' num2str(QuantileArray(i).EstimationQuantile*100,'%0.3f') ' %% Quantile:\t' num2str(QuantileArray(i).MassEstimate,'%0.2f') ' kg \n']);
    end
end
%% Calculate Collision Consequence
if isfield(AvailableServices,'CollisionConsequence') && AvailableServices.CollisionConsequence && ~isempty(PrimaryMass)
    % Prompt User to Input Secondary Mass (use Mass Estimation Quantiles if
    % available
    if ~exist('QuantileArray','var')
        % Set Default Dialog Box Text
        DialogText          = ['Input Secondary Object Mass Estimate For Collision Consequence Assessment - kg',...
                                ' (Default Corresponds to Overly Conservative Mass Estimate; ',...
                                'it is recommended that the user use the "Secondary Mass Estimation" module in',...
                                'conjunction with the "Collision Consequence" Module to better approximate secondary object mass)'];
        % Set Default Mass (DefaultQuantile)
        [MassQuantiles,MassQuantileArray] = ParseMassQuantileFlatFile(MassEstimateFile);
        try % Attempt to find Secondary Object in parsed file
            SecondaryIndex  = find(MassQuantileArray(:,1)==str2double(CDM.OBJECT2_OBJECT_DESIGNATOR));
            if isempty(SecondaryIndex) % If entry is not found, use a default estimate
                MassEstimate= 500;
            else
                % Interpolate mass estimate to desired quantile
                MassEstimate= interp1(MassQuantiles,MassQuantileArray(SecondaryIndex,2:end),DefaultQuantile,'pchip','extrap');
                DialogText  = ['Input Secondary Object Mass Estimate For Collision Consequence Assessment - kg (Default Corresponds to ' num2str(DefaultQuantile*100,'%0.2f') '% Quantile Mass Estimate)'];
            end
        catch
            MassEstimate    = 500;
        end
        % Request User Input Desired Secondary Mass
        answer          = inputdlg(DialogText...
                                    ,'Secondary Mass Estimate',1,{num2str(MassEstimate,'%0.2f')});
    else
        % Get Default Mass (Interpolate as needed)
        MassEstimate = interp1([QuantileArray(:).EstimationQuantile],[QuantileArray(:).MassEstimate],DefaultQuantile,'pchip','extrap');
        % Request User Input Desired Secondary Mass
        answer          = inputdlg(['Input Secondary Object Mass Estimate For Collision Consequence Assessment - kg (Default Corresponds to ' num2str(DefaultQuantile*100,'%0.2f') '% Quantile Mass Estimate)']...
                                    ,'Secondary Mass Estimate',1,{num2str(MassEstimate,'%0.2f')});
    end

    if isempty(answer)
        fprintf('\nNo Secondary Mass Estimate Input, Program will not Evaluate Collision Consequence\n')
        SecondaryMass         = answer;
    else
        SecondaryMass        = str2double(answer{1});
        if isnan(SecondaryMass)
            fprintf('\nNo Secondary Mass Estimate Input, or entry was ill-formatted, Program will not Evaluate Collision Consequence\n')
            SecondaryMass         = [];
        end
    end
    
    % Skip Processing if No Secondary Mass Input or No Primary Mass Input
    if ~(isempty(SecondaryMass) || isempty(PrimaryMass))
        % Print Spacer
        fprintf('\n');
        
        % Get Relative Velocity Between Objects
        VRel = norm(v1_J2K*1000-v2_J2K*1000);

        % Evaluate Collision Consequence
        [Catastrophic,NumOfPieces] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass);

        if Catastrophic
            fprintf(['Prospective Collision is expected to be:     "CATASTROPHIC" in Nature Using a Secondary Object Mass Estimate of: ' pad(num2str(SecondaryMass,'%0.2f'),8,'left') ' kg\n']);
        else
            fprintf(['Prospective Collision is expected to be: "NON-CATASTROPHIC" in Nature Using a Secondary Object Mass Estimate of: ' pad(num2str(SecondaryMass,'%0.2f'),8,'left') ' kg\n']);
        end
        % Set comand print line logic (digits if less than 1000, exponential if greater)
        if log10(NumOfPieces)<=3
            fprintf(['Expected Number of Debris Pieces if a Collision Occurs: ' pad(num2str(NumOfPieces,'%0.0f'),8,'left') '\n']);
        else
            fprintf(['Expected Number of Debris Pieces if a Collision Occurs: ' pad(num2str(NumOfPieces,'%0.2E'),8,'left') '\n']);
        end
    end
end
%% Calculate OD Quality Scores
if isfield(AvailableServices,'CDMSecondaryObjectODQA') && AvailableServices.CDMSecondaryObjectODQA
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

    % Calculate OD Quality Scores
    [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);

    % Collect Test Statistics
    ODQABinary = [CollectedScores(:,1); CollectedScores(~isnan(CollectedScores(:,2)),2); CollectedScores(~isnan(CollectedScores(:,3)),3); CollectedScores(~isnan(CollectedScores(:,4)),4)];
    temp = ODQAFlags(ODQABinary==0);

    % Print ODQA Results
    fprintf('\n');
    fprintf(['Orbit Determination Quality Score for Secondary Object "' num2str(CDM.OBJECT2_OBJECT_DESIGNATOR,'%05d') '":\t' num2str(CompositeScore,'%0.2f') '\n']);
    if ~isempty(temp)
        fprintf('The following Orbit Determination Quality Flags were triggered:\n');
        for i=1:length(temp)
            fprintf(['\t' temp{i} '\n'])
        end
    end
end
%% Monte Carlo Assessment of PC Using Equinoctial Sampling  
if (isfield(AvailableServices,'MonteCarloProbabilityofCollision') && AvailableServices.MonteCarloProbabilityofCollision)

    % Check if the 2D Probability of Collision has been calculated
    [Pc_Elrod,~,~,~]    = PcElrod(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
    % Check if the 3D Probability of Collision has been calculated
    [Pc_3D, Pc3D_Out]          = Pc3D_Hall(r1_J2K*1000,v1_J2K*1000,C1_J2K,r2_J2K*1000,v2_J2K*1000,C2_J2K,HBR);
    % Estimate Needed number of samples from 2DPc and 3DPc
    Nsample_kep = max([EstimateRequiredSamples(Pc_Elrod,Accuracy,MC_Confidence),...
                       EstimateRequiredSamples(Pc_3D,Accuracy,MC_Confidence)]);
%     Nsample_kep = 1E7;
    % Report Estimated Number of Trials Required
    fprintf(['\nDesired accuracy for Monte Carlo Simulation: ' num2str(Accuracy*100,'%0.1f') '%%\n']);
    fprintf(['Number of Monte Carlo samples required to estimate probability of collision to desired accuracy: \t' num2str(Nsample_kep,'%d') '\n'])

    % Bound Number of samples using minimum and maximum allowable sample counts
    Nsample_kep = max(1e2,Nsample_kep);
    Nsample_kep = min(1e8,Nsample_kep);
    % Report Number of trials to be performed
    fprintf(['Number of Monte Carlo samples to be performed: \t' num2str(Nsample_kep,'%d') '\n'])

    % Determine Batch Size min 1000, max 5000
    p = gcp;
    Nsample_batch = max([min([ceil(Nsample_kep/p.NumWorkers/1000)*1000 5000]) 1000]);

    % Get Primary Objecty Orbital Period
    a = -GM/2/(norm(v1_J2K*1000)^2/2-GM/norm(r1_J2K*1000));
    OrbitalPeriod = 2*pi()*sqrt(a^3/GM);

    % Get Bounds
    [tau0,tau1] = conj_bounds_Coppola(1e-16, HBR, (r2_J2K'-r1_J2K')*1000, (v2_J2K'-v1_J2K')*1000, C2_J2K+C1_J2K);
    tau0 = min([Pc3D_Out.Tmin tau0]);
    tau1 = max([Pc3D_Out.Tmax tau1]);
%     tau0 = Pc3D_Out.Tmin;
%     tau1 = Pc3D_Out.Tmax;
    % Report Bounds
    fprintf(['Coppola conjunction duration bounds (s): Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])

    % Check if event has low relative velocity
%     if (tau1-tau0) > LowVelThreshold*OrbitalPeriod
%         warning(['Event is classified as a low speed event: Coppola encounter time exceeds ' num2str(LowVelThreshold*100,'%.1f') '% of the Primary Object''s orbital period']);
%         warning(['Input encounter time will be modified to be no more than ' num2str(LowVelThreshold*100,'%.1f') '% of the orbital period'])
%         tau1    = LowVelThreshold*OrbitalPeriod*tau1/(tau1-tau0);
%         tau0    = tau1-LowVelThreshold*OrbitalPeriod;
%         fprintf(['Conjunction duration bounds (s) modified to: Dur = ' num2str(tau1-tau0,'%0.4f') ' (' num2str(tau0,'%0.4f') ' < t-TCA < ' num2str(tau1,'%0.4f') ')\n\n'])
%     end
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
    % E1TCA   = Cart2Eq(X1TCA); % Equinoctial elements at TCA
    % Jctoe1  = J_CART2EQ_Analytic(X1TCA); % Jacobian going from cartesian to equinoctial

    [KEP] = Cart2Kep([r1_J2K v1_J2K],'Mean','Rad'); % Get Keplerian elements to determine need for retrograde factor
    if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
        fr      = -1;
    else
        fr      = 1;
    end
    [~,n,af,ag,chi,psi,lM,~] = convert_cartesian_to_equinoctial(r1_J2K,v1_J2K,fr); % Equinoctial elements at TCA
    E1TCA   = [n af ag chi psi lM fr];
    Jctoe1 = jacobian_equinoctial_to_cartesian(E1TCA,X1TCA,fr); % Jacobian going from equinoctial to cartesian
    Jctoe1 = inv(Jctoe1); % Jacobian going from cartesian to equinoctial

    PEq1TCA = Jctoe1 * (C1_J2K/1e6) * Jctoe1'; % Equinoctial covariance at TCA

    X2TCA   = [r2_J2K v2_J2K]; % Cartesian state at TCA (km)
    % E2TCA   = Cart2Eq(X2TCA); % Equinoctial elements at TCA
    % Jctoe2  = J_CART2EQ_Analytic(X2TCA); % Jacobian going from cartesian to equinoctial

    [KEP] = Cart2Kep([r2_J2K v2_J2K],'Mean','Rad'); % Get Keplerian elements to determine need for retrograde factor
    if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
        fr      = -1;
    else
        fr      = 1;
    end
    [~,n,af,ag,chi,psi,lM,~] = convert_cartesian_to_equinoctial(r2_J2K,v2_J2K,fr); % Equinoctial elements at TCA
    E2TCA   = [n af ag chi psi lM fr];
    Jctoe2 = jacobian_equinoctial_to_cartesian(E2TCA,X2TCA,fr); % Jacobian going from equinoctial to cartesian
    Jctoe2 = inv(Jctoe2); % Jacobian going from cartesian to equinoctial

    PEq2TCA = Jctoe2 * (C2_J2K/1e6) * Jctoe2'; % Equinoctial covariance at TCA

    % Equinoctial Inputs
    [Pc_EQN_MC_all, Uc_EQN_MC_all, Nc_EQN_MCtfc_all,...
         Pc_EQN_MC_att, Uc_EQN_MC_att, Nc_EQN_MCtfc_att, ...
         Pc_EQN_MC_0, Uc_EQN_MC_0, Nc_EQN_MCtfc_0] = Pc_MC_Kep2body_parallel(...
            Nsample_kep, tfc_bins,...
            E1TCA', [], PEq1TCA,...
            E2TCA', [], PEq2TCA,...
            HBR, GM, 'k2bpri',...
            MC_Confidence, Nsample_batch,...
            FiguresFolder);


%     % Rectilinear Inputs
%     [~, ~, ~,...
%          Pc_CART_MC_att, Uc_CART_MC_att, Nc_CART_MCtfc_att,...
%          ~  , ~  , ~] =  Pc_MC_Kep2body_parallel(...
%             Nsample_kep, tfc_bins,...
%             X1TCA(1:3)'*1000, X1TCA(4:6)'*1000, C1_J2K,...
%             X2TCA(1:3)'*1000, X2TCA(4:6)'*1000, C2_J2K,...
%             HBR, GM, 'k2bpri',...
%             MC_Confidence, Nsample_batch,...
%             FiguresFolder);

    % Report Output Monte Carlo Statistics (EQN Sampling)
    fprintf('\n2-Body Propagation to Determine Probability of Collision - Equinoctial Inputs\n');
    fprintf(['\tNumber of Trials Resulting in Collision: \t' num2str(sum(Nc_EQN_MCtfc_att),'%d') '\n']);
    fprintf(['\tMonte Carlo Estimate of Collision Probability: \t' num2str(Pc_EQN_MC_att,'%0.2E') '\n']);
    fprintf(['\tConfidence Interval with ' num2str(MC_Confidence*100,'%0.1f') '%% Certainty: \t' num2str(Uc_EQN_MC_att(1),'%0.2E') ' < MC_Pc < ' num2str(Uc_EQN_MC_att(2),'%0.2E') '\n']);
    % Report Output Monte Carlo Statistics (Cartesian Sampling)
%     fprintf('\n2-Body Propagation to Determine Probability of Collision - Cartesian Inputs\n');
%     fprintf(['\tNumber of Trials Resulting in Collision: \t' num2str(sum(Nc_CART_MCtfc_att),'%d') '\n']);
%     fprintf(['\tMonte Carlo Estimate of Collision Probability: \t' num2str(Pc_CART_MC_att,'%0.2E') '\n']);
%     fprintf(['\tConfidence Interval with ' num2str(MC_Confidence*100,'%0.1f') '%% Certainty: \t' num2str(Uc_CART_MC_att(1),'%0.2E') ' < MC_Pc < ' num2str(Uc_CART_MC_att(2),'%0.2E') '\n']);
end
    %% Figure Generation for Analysis Purposes (Specifically recurring events)
%     [Pcdotmd_kep,Pcdotlo_kep,Pcdothi_kep] = ...
%         MC_Pc_limits(Nc_EQN_MCtfc_att,Nsample_kep,-MC_Confidence);
%     
%     % Calculate the cummulative MC Nc and Pc values for each bin
%     
%     Pc_kep_tfc_bins = Nc_EQN_MCtfc_att/Nsample_kep;
% 
%     Pcdot_kep_tfc_bins = Pcdotmd_kep/tfc_bins.del;
%     Pcdlo_kep_tfc_bins = Pcdotlo_kep/tfc_bins.del;
%     Pcdhi_kep_tfc_bins = Pcdothi_kep/tfc_bins.del;
% 
%     Pc_kep_tfc_bins_cum = cumsum(Pc_kep_tfc_bins);
%     
%     plot(tfc_bins.xmd,Pcdot_kep_tfc_bins, ...
%         'Marker','o', ...
%         'MarkerFaceColor','k', ...
%         'MarkerEdgeColor','k', ...
%         'MarkerSize',5, ...
%         'LineStyle','none', ...
%         'Color','k');
%     hold on;
%     ymn = Inf; ymx = -Inf;
%     for n=1:tfc_bins.N
%         xplt = [tfc_bins.xmd(n) tfc_bins.xmd(n)];
%         yplt = [Pcdlo_kep_tfc_bins(n) Pcdhi_kep_tfc_bins(n)];
%         plot(xplt,yplt,'-k');
%         ymn = min(ymn,yplt); ymx = max(ymx,yplt);
%     end
%     
%     if exist('Nc_CART_MCtfc_att','var')
%         [Pcdotmd_kep,Pcdotlo_kep,Pcdothi_kep] = ...
%             MC_Pc_limits(Nc_CART_MCtfc_att,Nsample_kep,-MC_Confidence);
% 
%         % Calculate the cummulative MC Nc and Pc values for each bin
% 
%         Pc_kep_tfc_bins = Nc_EQN_MCtfc_att/Nsample_kep;
% 
%         Pcdot_kep_tfc_bins = Pcdotmd_kep/tfc_bins.del;
%         Pcdlo_kep_tfc_bins = Pcdotlo_kep/tfc_bins.del;
%         Pcdhi_kep_tfc_bins = Pcdothi_kep/tfc_bins.del;
% 
%         Pc_kep_tfc_bins_cum = cumsum(Pc_kep_tfc_bins);
% 
%         plot(tfc_bins.xmd,Pcdot_kep_tfc_bins, ...
%             'Marker','s', ...
%             'MarkerFaceColor','r', ...
%             'MarkerEdgeColor','r', ...
%             'MarkerSize',5, ...
%             'LineStyle','none', ...
%             'Color','r');
%         hold on;
%         ymn = Inf; ymx = -Inf;
%         for n=1:tfc_bins.N
%             xplt = [tfc_bins.xmd(n) tfc_bins.xmd(n)];
%             yplt = [Pcdlo_kep_tfc_bins(n) Pcdhi_kep_tfc_bins(n)];
%             plot(xplt,yplt,'-r');
%             ymn = min(ymn,yplt); ymx = max(ymx,yplt);
%         end
%     end
%     hold off;
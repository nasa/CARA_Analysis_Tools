function [Analysis,SigData,RegData,AllData] = BFMCMultiCAProcess(conjID,inputPath,outputPath0)
% BFMCMultiCAProcess - Use BFMC long duration relative distance minima
%                      table (i.e., *.min) to calculate Nc and Pc for a
%                      multi-encounter interaction.
%                      (For CARA analysis team internal use)
%
% Syntax: [Analysis, SigData, RegData, AllData] = BFMCMultiCAProcess(conjID, inputPath, outputPath0);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    conjID     -   Conjunction ID string
%
%    inputPath  -   Path to BFMC run output folder containing *.eci and
%                   *.min files
%
%    outputPath0 -  Path to output directory for results and plots
%
% =========================================================================
%
% Output:
%
%   Analysis   -   Table containing analysis recommendation
%
%   SigData    -   Table of significant approach events
%                  (Pc >= sigPcLevel)
%
%   RegData    -   Table of all approach events that registered a Pc result
%                  (Pc >= 0)
%
%   AllData    -   Table of all approach events
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------
    
    % Initializations and parameters
    
    % Initialize output
    Analysis = []; SigData = []; RegData = []; AllData = [];

    % Run parameters
    verbose = true; debug_plotting = false;
    
    % Plotting option: 0 = no plots
    %                  1 = only multi-conjunction plots
    %                  2 = also individual conjunction plots
    plotting_option = 2;
    
    % Mahalanobis distance plotting option: 1 => plot MD & real(sqrt(MS2)) 
    %                                     : 2 => plot MD2 & MS2
    %                                     (make negative for 3DNc MDeff)
    MD_plotting_option = -1;
    
    % Writing option: 0 = no csv files output
    %                 1 = only analysis recommendation csv file output
    %                 2 = also significant event csv file output
    %                 3 = also registered event csv file output
    %                 4 = also all events csv file output
    writing_option = 4;
    
    % Green yellow and red Pc levels
    gPcLevel = 1e-10;
    yPcLevel = 1e-7;
    rPcLevel = 1e-4;
    
    % Minimum Pc value that is considered significant
    sigPcLevel = gPcLevel;

    % Minimum Pc value that to show on multi-conjunction plot
    % pltPcLevel = gPcLevel;
    pltPcLevel = 1e-12;
    
    % Cutoff level for HBR-adjusted Mahalanobis distances (MDcut)
    MScut = 1e3; MS2cut = MScut^2;
    
    % Cutoff levels for NPD-affected and temporally-extended conjunctions
    NPDEffectCut = 0.05;    
    ExtendCut = 0.5;
    BlendCut = ExtendCut; LongCut = BlendCut;
    
    % TCA offset level that prompts an error (s)
    dTCAmax = 1e-3;

    % Earth gravitational constant (EGM-96) [km^3/s^2]
    % GM = 3.986004418e5;   
    
    % Set up for 3D-Pc calculations
    Pc3Dpar = default_params_Pc3D_Hall([],false);
    
    % Pc values <= this are considered zero    
    tinyPcLevel = Pc3Dpar.Pc_tiny; 
    
    % Eigenvalue clipping parameter
    Fclip = Pc3Dpar.Fclip;
    
    % Screening volume processing
    screening_volume_processing = true;
    
    % Misc
    Pcfmt = '%0.7e'; % Format for Pc displaying
    dsstring = 'yyyy-mm-dd HH:MM:SS.FFF'; % ISO 8601 date string format
    MDlabel = 'N{\sigma} Distance'; % Mahalanobis Distance
    MD2label = 'M^2'; 
    
    % Define the primary and seconary object number strings from the
    % conjunction ID string
    
    [priIDStr, rem] = strtok(conjID, '_');
    [cnjIDStr, rem] = strtok(rem,    '_');
    [secIDStr, rem] = strtok(rem,    '_');
    
    if ~strcmpi(cnjIDStr,'conj')
        warning('Invalid conjunction ID string');
        return;
    else
        if verbose
            disp(' '); disp(' ');
            disp(['Processing conjunction: ' conjID]);
        end
    end
    
    % Get the creation time string
    [TCAIDStA, rem] = strtok(rem,    '_'); %#ok<ASGLU>
    [TCAIDStB, rem] = strtok(rem,    '_'); %#ok<ASGLU>
    [TCRIDStA, rem] = strtok(rem,    '_');
    [TCRIDStB, ~  ] = strtok(rem,    '_');
    TCRIDStr = [TCRIDStA '_' TCRIDStB];
    
    % Ensure the output folder for this conjunction exists
    % outputPath = fullfile(outputPath0,conjID);
    outputPath = outputPath0;
    if exist(outputPath, 'dir') == 0
        mkdir(outputPath);
    end
    
    % Search for existing ECI ephemeris files
    
    [priIDEph,priEphemFile] = getEphemFile(priIDStr,inputPath);
    priEphemExists = ~isempty(priEphemFile);
    
    priID = str2double(priIDStr); if isnan(priID); priID = priIDStr; end
    secID = str2double(secIDStr); if isnan(secID); secID = secIDStr; end
    
    [secIDEph,secEphemFile] = getEphemFile(secIDStr,inputPath);
    secEphemExists = ~isempty(secEphemFile);
    
    % Get the VCM files used to create the ephemeris files
    
    [priVCMFile] = getVCMFile(priIDStr,inputPath);
    [secVCMFile] = getVCMFile(secIDStr,inputPath);

    if priEphemExists && secEphemExists
        if verbose
            disp('Primary and secondary ephemeris files found:');
            [~,fn1,ex1] = fileparts(priEphemFile);
            if isempty(priVCMFile)
                VCMstr = ' (VCM file not found)';
            else
                VCMstr = [' (' priVCMFile ')'];
            end
            disp(['  ' fn1 ex1 VCMstr]);
            if isempty(secVCMFile)
                VCMstr = ' (VCM file not found)';
            else
                VCMstr = [' (' secVCMFile ')'];
            end            
            [~,fn2,ex2] = fileparts(secEphemFile);
            disp(['  ' fn2 ex2 VCMstr]);
        end
    else
        warning('Cannot find primary and/or secondary ephemeris files');
        return;
    end
    
    % Get the HBR
   
    outputMatFile = fullfile(inputPath,'output.mat');
    HBRLoaded = false;
    if exist(outputMatFile,'file')
        variableInfo = who('-file',outputMatFile);
        if ismember('HBR',variableInfo)
            load(outputMatFile,'HBR');
            HBRLoaded = true;
        end
    end
    if ~HBRLoaded
        cshFile = fullfile(inputPath,'vcm_single_checkout.csh');
        fileLines = splitlines(fileread(cshFile));
        for i = 1:length(fileLines)
            if startsWith(fileLines(i),'set a13')
                HBRstr = strrep(fileLines(i),'set a13=''','');
                HBRstr = strrep(HBRstr,'''','');
                HBR = str2double(HBRstr);
                HBRLoaded = true;
                break;
            end
        end
    end
    if ~HBRLoaded
        error(['Could not find HBR for ' conjID]);
    end

    % HBR in m & km
    HBR_m = HBR;
    HBR = HBR / 1000;
    
    % Screening volume processing
    SVprocessing = screening_volume_processing;
    if SVprocessing
        
        % OLD METHOD
        % % Read screening volume table, if available
        % if ~exist('SVtable.mat','file')
        %     SVprocessing = false;
        %     warning('Screening volume table not found');
        % else
        %     SVtable = []; load('SVtable.mat');
        %     [SVfound,SVindex] = ismember(priID,SVtable.SatID);
        %     if ~SVfound
        %         SVprocessing = false;
        %         warning('Primary not found in screening volume table');
        %     else
        %         SVRICdims = [SVtable.SV_R_km(SVindex) ...
        %                      SVtable.SV_I_km(SVindex) ...
        %                      SVtable.SV_C_km(SVindex)];
        %         if any(isnan(SVRICdims)) || any(SVRICdims <= 0)
        %             SVprocessing = false;
        %             warning('Primary has invalid entries in screening volume table');
        %         end
        %     end
        % end
        
        % Get mission table
        GMTpars.GetScreeningVolumes = true;
        MissionTable = GetMissionTable('',GMTpars);
        if priID == 36131
            warning('Handling special case of GEO 36131 -- mapping to 27566');
            [SVfound,SVindex] = ismember(27566,MissionTable.ObjectID);
        else
            [SVfound,SVindex] = ismember(priID,MissionTable.ObjectID);
        end
        if ~SVfound
            SVprocessing = false;
            warning('Primary not found in screening volume table');
        else
            SVRICdims = [MissionTable.TaskingR(SVindex) ...
                         MissionTable.TaskingI(SVindex) ...
                         MissionTable.TaskingC(SVindex)];
            SVRICshape = MissionTable.TaskingShape{SVindex};
            if strcmpi(SVRICshape,'Sphere')
                SVRICdims(2:3) = SVRICdims(1);
            end
            if any(isnan(SVRICdims))             || ...
               any(SVRICdims <= 0)               || ...
               isempty(SVRICshape)               || ...
               ~(strcmpi(SVRICshape,'Ellipsoid') || ...
                 strcmpi(SVRICshape,'Box')       || ...
                 strcmpi(SVRICshape,'Sphere'))
                SVprocessing = false;
                warning('Primary has invalid entries in screening volume table');
            end
        end
        
    end
    
     % Eigenvalue clipping parameter (km^2)
    Lclip = (Fclip*HBR)^2;
        
    % Read the ephemeris tables
    
    if verbose
        disp('Reading standard ephemeris');
    end    
    
    varNames = {'epoch','rx','ry','rz','vx','vy','vz',...
        'c11','c21','c22','c31','c32','c33','c41','c42','c43','c44',...
        'c51','c52','c53','c54','c55','c61','c62','c63','c64','c65',...
        'c66'};

    formatSpec = ['%f ' repmat('%e ',[1 26]) '%e\n'];
    numElems = 28;
    
    priEphem = readFortranAsTable(priEphemFile,numElems,formatSpec,varNames);
    secEphem = readFortranAsTable(secEphemFile,numElems,formatSpec,varNames);
    
    % Extract the ephemeris times, and ensure they are the same in both
    ephTimes = priEphem.epoch;
    if ~isequal(ephTimes,secEphem.epoch)
        error('Epoch values do not match for all rows in ephemeris files!')
    end
    NephTimes = numel(ephTimes);
    
    % Pri and sec ephemeris times are measured in days after the reference
    % time of 1969-12-31 00:00:00.0
    
    dnReference = datenum([1969 12 31 0 0 0]);
    
    dnPropBeg = ephTimes(1)+dnReference;
    [dsFullBeg,dsPropBeg] = convert_dn_to_ds(dnPropBeg);
    
    dnPropEnd = ephTimes(end)+dnReference;
    [dsFullEnd,dsPropEnd] = convert_dn_to_ds(dnPropEnd);
    
    if verbose
        disp('Ephemeris parameters;');
        disp([' Span: ' dsFullBeg ' to ' dsFullEnd]);
        difTimes = diff(ephTimes);
        disp([' Min/med/max time steps (s): ' ...
              num2str(min(difTimes)*86400) '  ' ...
              num2str(median(difTimes)*86400) '  ' ...
              num2str(max(difTimes)*86400)]);
    end
    
    % Generate the output file root
    outputRoot = [priIDStr '_' secIDStr '_' dsPropBeg '_to_' dsPropEnd];
    
    % Generate the ECI ephemeris structure
    
    STDEPH.T  = ephTimes';
    STDEPH.X1 = NaN(6,NephTimes);
    STDEPH.P1 = NaN(6,6,NephTimes);
    STDEPH.X2 = NaN(6,NephTimes);
    STDEPH.P2 = NaN(6,6,NephTimes);
    
    for i = 1:NephTimes
        
        % ECI pos/vel vectors
        r1 = [priEphem.rx(i) priEphem.ry(i) priEphem.rz(i)];
        v1 = [priEphem.vx(i) priEphem.vy(i) priEphem.vz(i)];
        r2 = [secEphem.rx(i) secEphem.ry(i) secEphem.rz(i)];
        v2 = [secEphem.vx(i) secEphem.vy(i) secEphem.vz(i)];
        
        % RIC covariances
        C1 = [priEphem.c11(i) priEphem.c21(i) priEphem.c31(i) priEphem.c41(i) priEphem.c51(i) priEphem.c61(i);
              priEphem.c21(i) priEphem.c22(i) priEphem.c32(i) priEphem.c42(i) priEphem.c52(i) priEphem.c62(i);
              priEphem.c31(i) priEphem.c32(i) priEphem.c33(i) priEphem.c43(i) priEphem.c53(i) priEphem.c63(i);
              priEphem.c41(i) priEphem.c42(i) priEphem.c43(i) priEphem.c44(i) priEphem.c54(i) priEphem.c64(i);
              priEphem.c51(i) priEphem.c52(i) priEphem.c53(i) priEphem.c54(i) priEphem.c55(i) priEphem.c65(i);
              priEphem.c61(i) priEphem.c62(i) priEphem.c63(i) priEphem.c64(i) priEphem.c65(i) priEphem.c66(i)];
        C2 = [secEphem.c11(i) secEphem.c21(i) secEphem.c31(i) secEphem.c41(i) secEphem.c51(i) secEphem.c61(i);
              secEphem.c21(i) secEphem.c22(i) secEphem.c32(i) secEphem.c42(i) secEphem.c52(i) secEphem.c62(i);
              secEphem.c31(i) secEphem.c32(i) secEphem.c33(i) secEphem.c43(i) secEphem.c53(i) secEphem.c63(i);
              secEphem.c41(i) secEphem.c42(i) secEphem.c43(i) secEphem.c44(i) secEphem.c54(i) secEphem.c64(i);
              secEphem.c51(i) secEphem.c52(i) secEphem.c53(i) secEphem.c54(i) secEphem.c55(i) secEphem.c65(i);
              secEphem.c61(i) secEphem.c62(i) secEphem.c63(i) secEphem.c64(i) secEphem.c65(i) secEphem.c66(i)];
    
        STDEPH.X1(:,i) = [r1 v1]';
        STDEPH.P1(:,:,i) = RIC2ECI(C1,r1,v1);
        STDEPH.X2(:,i) = [r2 v2]';
        STDEPH.P2(:,:,i) = RIC2ECI(C2,r2,v2);

    end
    
    % Calculate the separation distances and the Mahalanobis distances for
    % the combined ephemeris
    if verbose
        disp('Calculating HBR-adjusted Mahalanobis distances for standard ephemeris');
    end
    [STDEPH.RD,STDEPH.MD2,STDEPH.MS2] = calc_sep_maha_distances(STDEPH,HBR,Lclip);
    
    % Number of points in standard ephemeris
    NstdTimes = numel(STDEPH.T);

    % Read the detected minima in rel.pos. magnitudes, each of which
    % represents an approach between to two satellites
    
    if verbose
        disp('Reading minima ephemeris');
    end
    
    varNames = {'epoch','rdiff','vdiff',...
        'r1x','r1y','r1z','v1x','v1y','v1z',...
        'c1_11','c1_21','c1_22','c1_31','c1_32','c1_33','c1_41','c1_42','c1_43','c1_44',...
        'c1_51','c1_52','c1_53','c1_54','c1_55','c1_61','c1_62','c1_63','c1_64','c1_65',...
        'c1_66',...
        'r2x','r2y','r2z','v2x','v2y','v2z',...
        'c2_11','c2_21','c2_22','c2_31','c2_32','c2_33','c2_41','c2_42','c2_43','c2_44',...
        'c2_51','c2_52','c2_53','c2_54','c2_55','c2_61','c2_62','c2_63','c2_64','c2_65',...
        'c2_66'};
    
    formatSpec = [repmat('%f ',[1 3]) repmat('%e ',[1 53]) '%e\n'];
    numElems = 57;
    
    minEphem = readFortranAsTable( ...
        fullfile(inputPath,[priIDEph '_' secIDEph '.min']), ...
        numElems,formatSpec,varNames);
    
    % Eliminate the end points in the minima ephemeris if
    % they coincide with the end points of the main ephemeris
    
    if height(minEphem) > 1 && minEphem.epoch(1) == ephTimes(1)
        minEphem = minEphem(2:end,:);
    end
    
    if height(minEphem) > 1 && minEphem.epoch(end) == ephTimes(end)
        minEphem = minEphem(1:end-1,:);
    end
    
    minTimes = minEphem.epoch;
    NminTimes = numel(minTimes);
    
    MINEPH.T  = minTimes';
    MINEPH.X1 = NaN(6,NminTimes);
    MINEPH.P1 = NaN(6,6,NminTimes);
    MINEPH.X2 = NaN(6,NminTimes);
    MINEPH.P2 = NaN(6,6,NminTimes);
    
    for i = 1:NminTimes

        % ECI pos/vel vectors
        r1 = [minEphem.r1x(i) minEphem.r1y(i) minEphem.r1z(i)];
        v1 = [minEphem.v1x(i) minEphem.v1y(i) minEphem.v1z(i)];
        r2 = [minEphem.r2x(i) minEphem.r2y(i) minEphem.r2z(i)];
        v2 = [minEphem.v2x(i) minEphem.v2y(i) minEphem.v2z(i)];
        
        % RIC covariances
        
        C1 = [minEphem.c1_11(i) minEphem.c1_21(i) minEphem.c1_31(i) minEphem.c1_41(i) minEphem.c1_51(i) minEphem.c1_61(i);
              minEphem.c1_21(i) minEphem.c1_22(i) minEphem.c1_32(i) minEphem.c1_42(i) minEphem.c1_52(i) minEphem.c1_62(i);
              minEphem.c1_31(i) minEphem.c1_32(i) minEphem.c1_33(i) minEphem.c1_43(i) minEphem.c1_53(i) minEphem.c1_63(i);
              minEphem.c1_41(i) minEphem.c1_42(i) minEphem.c1_43(i) minEphem.c1_44(i) minEphem.c1_54(i) minEphem.c1_64(i);
              minEphem.c1_51(i) minEphem.c1_52(i) minEphem.c1_53(i) minEphem.c1_54(i) minEphem.c1_55(i) minEphem.c1_65(i);
              minEphem.c1_61(i) minEphem.c1_62(i) minEphem.c1_63(i) minEphem.c1_64(i) minEphem.c1_65(i) minEphem.c1_66(i)];
        C2 = [minEphem.c2_11(i) minEphem.c2_21(i) minEphem.c2_31(i) minEphem.c2_41(i) minEphem.c2_51(i) minEphem.c2_61(i);
              minEphem.c2_21(i) minEphem.c2_22(i) minEphem.c2_32(i) minEphem.c2_42(i) minEphem.c2_52(i) minEphem.c2_62(i);
              minEphem.c2_31(i) minEphem.c2_32(i) minEphem.c2_33(i) minEphem.c2_43(i) minEphem.c2_53(i) minEphem.c2_63(i);
              minEphem.c2_41(i) minEphem.c2_42(i) minEphem.c2_43(i) minEphem.c2_44(i) minEphem.c2_54(i) minEphem.c2_64(i);
              minEphem.c2_51(i) minEphem.c2_52(i) minEphem.c2_53(i) minEphem.c2_54(i) minEphem.c2_55(i) minEphem.c2_65(i);
              minEphem.c2_61(i) minEphem.c2_62(i) minEphem.c2_63(i) minEphem.c2_64(i) minEphem.c2_65(i) minEphem.c2_66(i)];
        
        MINEPH.X1(:,i) = [r1 v1]';
        MINEPH.P1(:,:,i) = RIC2ECI(C1,r1,v1);
        MINEPH.X2(:,i) = [r2 v2]';
        MINEPH.P2(:,:,i) = RIC2ECI(C2,r2,v2);
        
    end
    
    % Calculate separation distances at minima
    if verbose
        disp('Calculating HBR-adjusted Mahalanobis distances for minima ephemeris');
    end
    [MINEPH.RD,MINEPH.MD2,MINEPH.MS2] = calc_sep_maha_distances(MINEPH,HBR,Lclip);
    
    if verbose
        disp('Approach or relative distance minima ephemeris parameters:');
        disp([' Number of minima: ' num2str(NminTimes)]);
        disp([' Minimum approach miss-distance = ' num2str(min(MINEPH.RD)) ' km']);
        disp([' Maximum approach miss-distance = ' num2str(max(MINEPH.RD)) ' km']);
    end
    
    clear priEphem secEphem minEphem minTimes ephTimes;
    
    % Create the combined ephemeris, with the standard and CA-minima
    % quantities combined and sorted in time.
    
    % Exclude ephemeris points that nearly coincide with the minima
    % times. This prevents 5-point state interpolation from failing due
    % to combined eph. times that are too close to one another
    StepFracMin = 0.1;
    excl_eph_time = false(size(STDEPH.T)); Nexcl = 0;
    for i = 1:NminTimes
        [c,b,coincidence] = find_bracketing_indices(MINEPH.T(i),STDEPH.T);
        if coincidence
            % Minimum perfectly coincides with an ephemeris point
            excl_eph_time(c) = true; Nexcl = Nexcl+1;
        else
            StepFrac = (MINEPH.T(i)-STDEPH.T(c))/(STDEPH.T(b)-STDEPH.T(c));
            if StepFrac < 0 || StepFrac > 1
                error('Invalid time-step fraction calculated');
            end
            if StepFrac < StepFracMin
                excl_eph_time(c) = true;  Nexcl = Nexcl+1;
            end
        end       
    end
    if verbose
        disp(['Number of eph times excluded as too close to minima: ' ...
            num2str(Nexcl)]);
    end
    % Combine the ephemeris tables, excluding eph points as required
    ndx = ~excl_eph_time;    
    CMBEPH.T   = cat( 2 , MINEPH.T   , STDEPH.T(ndx)      );
    CMBEPH.X1  = cat( 2 , MINEPH.X1  , STDEPH.X1(:,ndx)   );
    CMBEPH.P1  = cat( 3 , MINEPH.P1  , STDEPH.P1(:,:,ndx) );
    CMBEPH.X2  = cat( 2 , MINEPH.X2  , STDEPH.X2(:,ndx)   );
    CMBEPH.P2  = cat( 3 , MINEPH.P2  , STDEPH.P2(:,:,ndx) );
    CMBEPH.RD  = cat( 2 , MINEPH.RD  , STDEPH.RD(:,ndx)   );
    CMBEPH.MD2 = cat( 2 , MINEPH.MD2 , STDEPH.MD2(:,ndx)  );
    CMBEPH.MS2 = cat( 2 , MINEPH.MS2 , STDEPH.MS2(:,ndx)  );
    
    % Sort the combined ephemeris
    [CMBEPH.T,ndx] = sort(CMBEPH.T);
    CMBEPH.X1  = CMBEPH.X1(:,ndx);
    CMBEPH.P1  = CMBEPH.P1(:,:,ndx);
    CMBEPH.X2  = CMBEPH.X2(:,ndx);
    CMBEPH.P2  = CMBEPH.P2(:,:,ndx);
    CMBEPH.RD  = CMBEPH.RD(:,ndx);
    CMBEPH.MD2 = CMBEPH.MD2(:,ndx);
    CMBEPH.MS2 = CMBEPH.MS2(:,ndx);
    
    % Number of points in combined ephemeris
    NcmbTimes = numel(CMBEPH.T);
    
    % Define the relative ephemeris times (days)
    
    cmbTimesFirst = CMBEPH.T(1);
    % cmbTimesLast  = CMBEPH.T(end);

    cmbRelTimes = CMBEPH.T-cmbTimesFirst;
    minRelTimes = MINEPH.T-cmbTimesFirst;

    STDEPH.TREL = STDEPH.T-cmbTimesFirst;
    CMBEPH.TREL = CMBEPH.T-cmbTimesFirst;
    
    % Find the maxima of the relative position magnitudes to use as the
    % start and stop times of the sequence of conjunction segments
    
    if verbose
        disp('Finding maxima of separation distance curve');
    end

    [maxRelPosMags,imax,~,imin] = extrema(CMBEPH.RD,false,false); %#ok<ASGLU>   
    maxRelTimes = CMBEPH.T(imax)-cmbTimesFirst;
    
    % Add edge minima
    if excl_eph_time(1)
        imin = cat(2,1,imin);
    end
    if excl_eph_time(NstdTimes)
        imin = cat(2,imin,NstdTimes);
    end
    
    % Ensure that the minima were recovered for consistency
    if ~isequal(CMBEPH.RD(imin),MINEPH.RD) || ~isequal(CMBEPH.T(imin),MINEPH.T)
        error('Mismatch in determined relative distance minima');
    end
    
    % Define the ephemeris to be used for interpolation
    % INTEPH = STDEPH;
    INTEPH = CMBEPH;
    
    % Refine the minima of the M2HBR curve
    if verbose
        disp('Finding minima of HBR-adjusted Mahalanobis distance curve');
    end
    % % Bisect to at least 1/10 of eph step
    % Ttol = max(1e-4/86400,min(diff(STDEPH.T))/10);
    % Ttol = [Ttol NaN]; % Impose only absolute time tolerance
    Ttol = [NaN NaN]; % Impose no time tolerance
    MD2tol = [3e-2 1e-2]; % Impose both absolute and relative MD2 tolerance
    fun = @(T) interp_MS2(T,INTEPH,HBR,Lclip,0);
    extrema_types = 1; % Refine minima but not maxima
    endpoints = true; % Refine minima at endpoints
    rbeverbose = verbose;
    rbecheck = false;
    [TMS2min,MS2min,~,~,rbeconv,rbenbisect,TMD2buf,MS2buf,~,~] = ...
        refine_bounded_extrema(fun,INTEPH.TREL,INTEPH.MS2, ...
                               [],100,extrema_types,Ttol,MD2tol, ...
                               endpoints,rbeverbose,rbecheck);
    if ~rbeconv
        error(['MS2 minima search failed to converge in ' ...
            num2str(rbenbisect) ' bisections']);
    end
    
    % Interpolate MD2 not adjusted for HBR
    MD2buf = interp_MS2(TMD2buf,INTEPH,0,Lclip,0);

    % figure;
    % subplot(2,1,1)
    % plot(TMD2buf,sqrt(MS2buf),'-k');
    % hold on;
    % plot(cmbRelTimes,sqrt(CMBEPH.MS2),'.r');
    % plot(TMS2min,sqrt(MS2min),'ob');
    % hold off;
    % subplot(2,1,2);
    % plot(TMD2buf,sqrt(MD2buf),'-k');    
    % keyboard;
        
    % Set up for plotting    
    xfwt = 'bold';
    yfwt = 'bold';
    tfwt = 'bold';
    afwt = 'bold';
    alwd = 0.5;

    lnwdthin = 1;
    lnwdthick = lnwdthin+1;

    msiz = 4;
    
    gold = [1 0.667 0];
    drkgrn = [0 0.667 0];
    
    txtblack = '{0 0 0}';
    txtred = '{1 0 0}';
    txtblue = '{0 0 1}';
    txtgold = '{1 0.667 0}';
    
    gridmajcolor = [0.5 0.5 0.5];
    gridmincolor = [0.7 0.7 0.7];

    % Plot processing
    if plotting_option > 0
        figure;
        set(gcf,'Visible','off');
        p = get(gcf,'Position');
        % set(gcf,'Position',[p(1:2) 1075 800]);
        % movegui('center');
        % lft = p(1); bot = p(2); wid = p(3); hgt = p(4);        
        set(gcf,'position',[round(p(1:2)/4) round(p(3:4)*1.5)]);
        op = get(gcf,'OuterPosition');
        set(gcf,'InnerPosition',op);
        set(gcf,'Visible','on');
    end
    
    % Output data
    out.Pri = NaN(NminTimes,1);
    out.Sec = NaN(NminTimes,1);
    out.TCA = cell(NminTimes,1);
    out.HBR_m = NaN(NminTimes,1);
    out.MissDist = NaN(NminTimes,1);
    out.MissVrel = NaN(NminTimes,1);
    out.MD = NaN(NminTimes,1);
    out.MS = NaN(NminTimes,1);
    out.NminMS2 = NaN(NminTimes,1);
    out.Pc3D = NaN(NminTimes,1);
    out.Tpeak3D  = NaN(NminTimes,1);
    out.Ta3D = NaN(NminTimes,1);
    out.Tb3D = NaN(NminTimes,1);
    out.NPDEffect3D = NaN(NminTimes,1);
    out.Extend3D = NaN(NminTimes,1);
    out.Blend3D = NaN(NminTimes,1);
    out.Long3D = NaN(NminTimes,1);
    out.NmaxNcdot = NaN(NminTimes,1);    
    out.Pc2D = NaN(NminTimes,1);
    out.Tpeak2D  = NaN(NminTimes,1);
    out.Ta2D = NaN(NminTimes,1);
    out.Tb2D = NaN(NminTimes,1);
    out.Extend2D = NaN(NminTimes,1);
    out.Blend2D = NaN(NminTimes,1);
    out.Long2D = NaN(NminTimes,1);
    out.TAdatestr = cell(NminTimes,1);
    out.TBdatestr = cell(NminTimes,1);
    out.Conjunction_ID = cell(NminTimes,1);
    
    out.MissR = NaN(NminTimes,1);
    out.MissI = NaN(NminTimes,1);
    out.MissC = NaN(NminTimes,1);
    if SVprocessing
        out.InScr = NaN(NminTimes,1);
    end
    
    relTA = NaN(1,NminTimes);
    relTB = NaN(1,NminTimes);
    
    reltaua2D = NaN(1,NminTimes);
    reltaub2D = NaN(1,NminTimes);
    reltaua3D = NaN(1,NminTimes);
    reltaub3D = NaN(1,NminTimes);
    
    Tcum = []; Ncum = [];
    
    for i = 1:NminTimes

        % Extract TCA conjunction data
        r1 = MINEPH.X1(1:3,i)';
        v1 = MINEPH.X1(4:6,i)';
        r2 = MINEPH.X2(1:3,i)';
        v2 = MINEPH.X2(4:6,i)';
        C1ECI = MINEPH.P1(:,:,i);
        C2ECI = MINEPH.P2(:,:,i);

        % Calculate RIC miss components
        % Unit vectors in the primary's RIC directions
        h1   = cross(r1,v1);
        rhat = r1 / norm(r1);
        chat = h1 / norm(h1);
        ihat = cross(chat,rhat);
        % RIC to ECI rotation matrix
        ECItoRIC = [rhat; ihat; chat];
        % Primary-to-Secondary miss vector in primary's RIC frame
        RICmiss = ECItoRIC * (r2-r1)';
        % Output RIC miss components
        out.MissR(i) = RICmiss(1);
        out.MissI(i) = RICmiss(2);
        out.MissC(i) = RICmiss(3);

        % Screening volume processing
        if SVprocessing
            SVshape = SVRICshape;
            SVdims  = SVRICdims;
            if strcmpi(SVRICshape,'Sphere')
                SVshape = 'Ellipsoid';
                SVdims(:) = SVdims(1);
            end
            if strcmpi(SVshape,'Box')
                out.InScr(i) = abs(RICmiss(1)) <= SVdims(1) & ...
                               abs(RICmiss(2)) <= SVdims(2) & ...
                               abs(RICmiss(3)) <= SVdims(3);
            elseif strcmpi(SVshape,'Ellipsoid')
                rr = RICmiss./reshape(SVdims,size(RICmiss));
                out.InScr(i) = sum(rr.^2) < 1;
            else
                error('Invalid SVshape');
            end
        end

        % Set up debug plotting
        if debug_plotting
            clf; %#ok<UNRCH>
            subplot(3,1,1);
            plot(cmbRelTimes,CMBEPH.RD,'k');
            hold on;
            xmn = minRelTimes(i);
            xmx = minRelTimes(i);
            plot(minRelTimes(i),MINEPH.RD(i),'vr');
        end
        
        % Find the bounding maxima times for the conjunction segment
        % defined by this min-rel-distance approach
        iA = find(maxRelTimes < minRelTimes(i));
        if isempty(iA)
            if i == 1
                TA = cmbRelTimes(1);
            else
                error('Failed to find begin time of conjunction segment');
            end
        else
            iA = iA(end);
            TA = maxRelTimes(iA);
            if debug_plotting
                xmn = TA; %#ok<UNRCH>
                plot(TA,maxRelPosMags(iA),'^b');
            end
        end
        iB = find(maxRelTimes > minRelTimes(i));
        if isempty(iB)
            if i == NminTimes
                TB = cmbRelTimes(NcmbTimes);
            else
                error('Failed to find end time of conjunction segment');
            end
        else
            iB = iB(1);
            TB = maxRelTimes(iB);
            if debug_plotting
                xmx = TB; %#ok<UNRCH>
                plot(TB,maxRelPosMags(iB),'^b');
                hold off;
            end
        end
        
        % Find the number of MS2 minima within this conjunction segment
        ndx = (TA <= TMS2min) & (TMS2min <= TB) & ...
              (MS2min <= MS2cut);
        NminMS2 = sum(ndx);
        
        % Find the Mahalonobis distances within this conjunction segment
        NMD2buf = numel(TMD2buf);
        ndx = (TA <= TMD2buf) & (TMD2buf <= TB);
        MD = sqrt(max(0,min(MD2buf(ndx))));
        MS = sqrt(max(0,min(MS2buf(ndx))));
        
        if verbose
            disp(['CA ' num2str(i) ...
                '  MD = ' num2str(MD) '   MS = '  num2str(MS)]);
        end

        % Store output data
        out.Pri(i) = priID;
        out.Sec(i) = secID;
        dnTCA = MINEPH.T(i)+dnReference;        
        out.TCA{i} = [' ' datestr(dnTCA,dsstring)];
        % out.TCA{i} = datestr(dnTCA,dsstring);
        out.HBR_m(i) = HBR_m;
        out.MissDist(i) = MINEPH.RD(i);
        out.MissVrel(i) = norm(MINEPH.X2(4:6,i)-MINEPH.X1(4:6,i));
        out.MD(i) = MD;
        out.MS(i) = MS;
        out.NminMS2(i) = NminMS2;
        out.Conjunction_ID{i} = conjID;
        
        % Process the conjunctions that make the MD cut
        
        if MS <= MScut
            
            % Process approaches that make the MD-based cutoff
            
            if debug_plotting
                xrng = plot_range([xmn,xmx],0.10); %#ok<UNRCH>
                xlim(xrng);
                ymn = min([MS2buf(ndx) MD2buf(ndx)]);
                ymx = max([MS2buf(ndx) MD2buf(ndx)]);
                ymx = min(ymx,MScut);
                yrng = plot_range([ymn ymx],0.10);
                subplot(3,1,2)
                plot(TMD2buf,MS2buf,':r');
                hold on;
                plot(TMD2buf(ndx),MS2buf(ndx),'-r');
                plot(TMD2buf,MD2buf,':b');
                plot(TMD2buf(ndx),MD2buf(ndx),'--b');
                plot(xrng,[0 0],'-k');
                plot(xrng,[MScut MScut],'--k');
                hold off;
                xlim(xrng); ylim(yrng);
            end

            % Pc2D calcuations use the corrected-for-TCA-offsets states
            [dTCA,X1CA,X2CA] = FindNearbyCA([r1 v1]',[r2 v2]'); 
            if (dTCA > dTCAmax)
                error('TCA offset exceeds max alowable value');
            end
            [Pc2D,~,~,~] = Pc2D_Foster(X1CA(1:3)',X1CA(4:6)',C1ECI, ...
                                       X2CA(1:3)',X2CA(4:6)',C2ECI, ...
                                       HBR,1e-8,'circle');
            if verbose
                disp([' 2D-Pc = ' num2str(Pc2D,Pcfmt)]);
            end
            
            % Find the Coppola linear conjunction bounds
            r = X2CA(1:3)-X1CA(1:3);
            v = X2CA(4:6)-X1CA(4:6);
            C = C2ECI+C1ECI;
            [tau0,tau1] = conj_bounds_Coppola([1e-6 1e-16], HBR, r, v, C);
            if (Pc2D <= tinyPcLevel); Pc2D = 0; end
            
            % Pc2D conjunction time limits are the Coppola bounds defined
            % at gamma = 1e-6, in order to be comparable to the
            % Pc3D time limits later estimated at the same gamma level
            out.Ta2D(i) = tau0(1);
            out.Tb2D(i) = tau1(1);
            out.Tpeak2D(i) = 0.5*(out.Ta2D(i)+out.Tb2D(i));

            % Conjunction segment time limits relative to TCA (s)
            Pc3Dpar.Tmin_limit = (TA-minRelTimes(i))*86400;
            Pc3Dpar.Tmax_limit = (TB-minRelTimes(i))*86400;
            
            % Find the time limits that MS is below the cutoff limit
            idx = find(ndx & (MS2buf <= MS2cut));
            minidx = max(idx(1)-1,1);
            maxidx = min(idx(end)+1,NMD2buf);
            mincut = TMD2buf(minidx);
            maxcut = TMD2buf(maxidx);
            if debug_plotting
                hold on; %#ok<UNRCH>
                plot([mincut mincut],yrng,'--c');
                plot([maxcut maxcut],yrng,'--c');
                hold off;
            end

            % Set initial conjunction duration limits for Pc3D
            if (minidx == 1) && (maxidx == NMD2buf)
                % Initial Nc3D time limits are full segment durations
                Pc3Dpar.Tmin_initial = -Inf;
                Pc3Dpar.Tmax_initial =  Inf;
            else
                % Initial conjunction duration for Pc3D
                taua = (mincut-minRelTimes(i))*86400;
                taub = (maxcut-minRelTimes(i))*86400;
                % Expand the empirically estimated Pc3D time bounds to the 
                % Coppola conjunction bounds defined at gamma = 1e-16, to
                % help ensure that they are adequately wide
                Pc3Dpar.Tmin_initial = max(min(taua,tau0(2)),Pc3Dpar.Tmin_limit);
                Pc3Dpar.Tmax_initial = min(max(taub,tau1(2)),Pc3Dpar.Tmax_limit);
            end
            
            % Pc3Dpar.Tmin_initial = [];
            % Pc3Dpar.Tmax_initial = [];
                
            % Pc3D calculations can use the original TCA states  
            % Pc3Dpar.debug_plotting = 1;
            % Pc3Dpar.debug_plotting = 2*(Pc2D > 1e-10);
            [Pc3D, Pc3Dout] = Pc3D_Hall(r1*1e3,v1*1e3,C1ECI*1e6, ...
                                        r2*1e3,v2*1e3,C2ECI*1e6, ...
                                        HBR*1e3,Pc3Dpar);
            if ~Pc3Dout.converged
                
                Pc3D = -1;
                if verbose
                    disp(' 3D-Pc not converged');
                end
                
            else
                
                if verbose
                    disp([' 3D-Pc = ' num2str(Pc3D,Pcfmt) ...
                          ' MDmin = ' num2str(sqrt(Pc3Dout.MD2min)) ...
                          ' Neph = ' num2str(Pc3Dout.Neph) ...
                          ' Nrefine = ' num2str(Pc3Dout.nrefine)]);
                end
                if (Pc3D <= tinyPcLevel); Pc3D = 0; end
                
                % Calculate remediated 3D-Pc estimate, if necessary. This
                % should need to be done explicitly relatively rarely.
                if Pc3Dout.Qmean10RemStat || Pc3Dout.Qmean20RemStat
                    Pc3Dpar.remediate_NPD_TCA_eq_covariances = true;
                    Pc3Draw = Pc3D;
                    % Pc3Doutraw = Pc3Dout;
                    [Pc3Drem, Pc3Dremout] = ...
                        Pc3D_Hall(r1*1e3,v1*1e3,C1ECI*1e6, ...
                                  r2*1e3,v2*1e3,C2ECI*1e6, ...
                                  HBR*1e3,Pc3Dpar);
                    Pc3Dpar.remediate_NPD_TCA_eq_covariances = false;
                    if Pc3Dremout.converged
                        % Pc3D = Pc3Drem;
                        % Pc3Dout = Pc3Dremout;
                        out.NPDEffect3D(i) = 2*abs(Pc3Draw-Pc3Drem)/(Pc3Draw+Pc3Drem);                    
                        if verbose
                            disp([' Remediation required: ' num2str(out.NPDEffect3D(i))]);
                        end                        
                    else
                        if verbose
                            disp([' Remediation required and failed to converge']);
                        end
                        out.NPDEffect3D(i) = 1;
                    end
                    
                else
                    % Pc3Draw = Pc3D;
                    out.NPDEffect3D(i) = 0;
                end
                
            end

            if debug_plotting && Pc3Dout.converged
                subplot(3,1,3);
                tt = minRelTimes(i)+Pc3Dout.Teph/86400;
                ndx = Pc3Dout.Ncdot > 0;
                if sum(ndx) > 1
                    semilogy(tt,Pc3Dout.Ncdot,'-r');
                else
                    plot(tt,Pc3Dout.Ncdot,'-r');
                end
                hold on;
                plot(tt,Pc3Dout.Ncdot_SmallHBR,':k');
                hold off;
                xlim(xrng);
                titl = ['Nc = ' num2str(Pc3Dout.Nc,Pcfmt) ...
                      '  TaFrac = ' num2str(Pc3Dout.TaFrac) ...
                      '  TbFrac = ' num2str(Pc3Dout.TbFrac)];
                title(titl);
                if ~isnan(Pc3Dout.TaFrac) % && Pc3Dout.Nc > 1e-10
                    keyboard;
                end
            end

            % Define the time bounds for the full conjunction segment
            tbeg = Pc3Dout.Tmin_limit/86400;
            tend = Pc3Dout.Tmax_limit/86400;
            out.TAdatestr{i} = [' ' datestr(tbeg + dnTCA,dsstring)];
            out.TBdatestr{i} = [' ' datestr(tend + dnTCA,dsstring)];
            % out.TAdatestr{i} = datestr(tbeg + dnTCA,dsstring);
            % out.TBdatestr{i} = datestr(tend + dnTCA,dsstring);
            relTA(i) = tbeg + minRelTimes(i);
            relTB(i) = tend + minRelTimes(i);

            % Define the time bounds for the 2D-Pc encounter, and the
            % blending and long-duration fractions
            out.Pc2D(i) = Pc2D;
            tbeg = out.Ta2D(i)/86400;
            tend = out.Tb2D(i)/86400;
            reltaua2D(i) = tbeg + minRelTimes(i);
            reltaub2D(i) = tend + minRelTimes(i);
            out.Blend2D(i) = max(out.Ta2D(i)/Pc3Dpar.Tmin_limit, ...
                                 out.Tb2D(i)/Pc3Dpar.Tmax_limit);
            out.Long2D(i) = (out.Tb2D(i)-out.Ta2D(i))/ ...
                (Pc3Dpar.Tmax_limit-Pc3Dpar.Tmin_limit);
            out.Extend2D(i) = max(out.Blend2D(i),out.Long2D(i));

            % Define the time bounds for the 3D-Pc encounter, and the
            % blending and long-duration fractions
            out.Pc3D(i) = Pc3D;
            out.Tpeak3D(i) = Pc3Dout.TpeakConj; 
            out.Ta3D(i) = Pc3Dout.TaConj;
            out.Tb3D(i) = Pc3Dout.TbConj;
            tbeg = out.Ta3D(i)/86400;
            tend = out.Tb3D(i)/86400;
            reltaua3D(i) = tbeg + minRelTimes(i);
            reltaub3D(i) = tend + minRelTimes(i);
            out.Blend3D(i) = max(out.Ta3D(i)/Pc3Dpar.Tmin_limit, ...
                                 out.Tb3D(i)/Pc3Dpar.Tmax_limit);
            out.Long3D(i) = (out.Tb3D(i)-out.Ta3D(i))/ ...
                (Pc3Dpar.Tmax_limit-Pc3Dpar.Tmin_limit);
            out.Extend3D(i) = max(out.Blend3D(i),out.Long3D(i));
            
            % Include all maxima here (rather than just significant maxima)
            % because 2D-Pc may fix onto the wrong single Ncdot maximum an
            % get inaccurate results
            out.NmaxNcdot(i) = Pc3Dout.Ncmaxima;

            % Calculate the cumulative Nc curve if the 3D-Pc algorithm
            % converged
            if Pc3D > 0
                if isempty(Ncum)
                    Tcum = minRelTimes(i)+calc_cum_times(Pc3Dout.Teph)/86400;
                    Ncum = Pc3Dout.Nccum;
                else
                    Tcum = cat(2,Tcum,minRelTimes(i)+calc_cum_times(Pc3Dout.Teph)/86400);
                    Ncum = cat(2,Ncum,Ncum(end)+Pc3Dout.Nccum);
                end
            end
            
            % Make the event plot, but only for cases where 3DPc converged
            % and has a Pc value sufficiently large
            
            if (plotting_option >= 1) && (Pc3D >= pltPcLevel)
                
                % Clear plot
                clf;
                
                % Title info
                subplot(32,2,[15 16]);
                % plot([NaN NaN],[NaN NaN]);
                axis off;
                
                titl = [];
                
                ttl = ' ';
                titl = cat(1,titl,{ttl});

                ttl = ['Primary = ' priIDStr ...
                    '   Secondary = ' secIDStr ...
                    '   HBR = ' num2str(HBR*1e3) 'm   (' ...
                    'approach ' num2str(i) ' of ' num2str(NminTimes) ')'];
                titl = cat(1,titl,{ttl});
                
                ttl = ' ';
                titl = cat(1,titl,{ttl});
                
                % % Write propagation time info
                % ttl = ['Ephemeris start: ' dsFullBeg];
                % titl = cat(1,titl,{ttl});
                % ttl = ['Ephemeris end:   ' dsFullEnd];
                % titl = cat(1,titl,{ttl});

                % Write closest approach info
                [TCAstr,TCAtrm] = convert_dn_to_ds(MINEPH.T(i)+dnReference);
                dca = MINEPH.RD(i);
                if dca < 1
                    dcastr = [num2str(dca*1e3,'%0.1f') ' m'];
                else
                    if dca < 10
                        dcastr = [num2str(dca,'%0.3f') ' km'];
                    else
                        dcastr = [num2str(dca,'%0.2f') ' km'];
                    end
                end
                ttl = ['Miss Distance = ' dcastr '   TCA @ ' TCAstr];
                titl = cat(1,titl,{ttl});
                
                ttl = ' ';
                titl = cat(1,titl,{ttl});

                % Write Pc info
                [~,ndx] = max(Pc3Dout.Ncdot);
                TPKstr = convert_dn_to_ds( ...
                    Pc3Dout.Teph(ndx)/86400+MINEPH.T(i)+dnReference);
                Pcstr = num2str(Pc3D,'%0.3e');
                ttl = ['Nc = ' Pcstr '   Peak Rate @ ' TPKstr];
                if (Pc3D < yPcLevel)
                    txtclr = '{0 0.667 0}';
                    Ncdclr = drkgrn;
                elseif (Pc3D < rPcLevel)
                    txtclr = '{1 0.667 0}';
                    Ncdclr = gold;
                else
                    txtclr = '{1 0 0}';
                    Ncdclr = 'r';
                end
                ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
                titl = cat(1,titl,{ttl});
                
                ttl = ' ';
                titl = cat(1,titl,{ttl});
                
                % Usage indicators
                
                ttlN = ['NPDEffect=' smart_exp_format(out.NPDEffect3D(i),2)];
                if out.NPDEffect3D(i) > NPDEffectCut
                    txtclr = txtred;
                else
                    txtclr = txtblue;                    
                end
                ttlN = ['{\color[rgb]' txtclr ttlN '}']; %#ok<AGROW>
                
                ttlL = ['Duration=' smart_exp_format(out.Long3D(i),2)];
                if out.Long3D(i) > LongCut
                    txtclr = txtred;
                else
                    txtclr = txtblue;                    
                end
                ttlL = ['{\color[rgb]' txtclr ttlL '}']; %#ok<AGROW>
                
                ttlB = ['Blending=' smart_exp_format(out.Blend3D(i),2)];
                if out.Blend3D(i) > BlendCut
                    txtclr = txtred;
                else
                    txtclr = txtblue;                    
                end
                ttlB = ['{\color[rgb]' txtclr ttlB '}']; %#ok<AGROW>
                
                ttlV = ['NumRatePeaks=' num2str(out.NmaxNcdot(i))];
                if out.NmaxNcdot(i) > 1
                    txtclr = txtgold;
                else
                    txtclr = txtblue;                    
                end
                ttlV = ['{\color[rgb]' txtclr ttlV '}']; %#ok<AGROW>
                
                ttl = ['Usage Indicators: ' ttlN ', ' ttlL ', ' ttlB ', ' ttlV];
                titl = cat(1,titl,{ttl});
                                 
                % ttl = ' ';
                % titl = cat(1,titl,{ttl});
                
                % Write assembled title information
                title(titl,'FontWeight',tfwt,'FontAngle','Italic');
                
                % Plot relative distance over the conjunction segment
                subplot(4,2,3);
                xmn = relTA(i);
                xmx = relTB(i);
                plot(cmbRelTimes,CMBEPH.RD,':k','LineWidth',lnwdthin);
                ndx_RD = (xmn <= cmbRelTimes) & (cmbRelTimes <= xmx);
                hold on;
                plot(cmbRelTimes(ndx_RD),CMBEPH.RD(ndx_RD),'-k','LineWidth',lnwdthick);
                % plot(minRelTimes(i),MINEPH.RD(i),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
                hold off;
                xrng = plot_range([xmn xmx],0.10);
                yrng = plot_range(CMBEPH.RD(ndx_RD),0.10,0.35);
                yrng(yrng < 0) = 0;
                xlim(xrng); ylim(yrng);
                % xlabl = 'Ephemeris Time (days)';
                xlabl = 'Time (days)';
                % xlabel(xlabl,'FontWeight',xfwt);
                ylabl = 'Separation (km)';
                ylabel(ylabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);
                
                % Plot the MD curves
                subplot(4,2,5);
                ndx_MD = (relTA(i) <= TMD2buf) & (TMD2buf <= relTB(i));
                xplt_MD = TMD2buf(ndx_MD);
                TMD23D = minRelTimes(i)+Pc3Dout.Teph/86400;
                ndx_3D = (relTA(i) <= TMD23D) & (TMD23D <= relTB(i));
                xplt_3D = TMD23D(ndx_3D);
                if abs(MD_plotting_option) == 1
                    yplt1_MD = real(sqrt(MD2buf(ndx_MD)));
                    yplt2_MD = real(sqrt(MS2buf(ndx_MD)));
                    yplt3_MD = real(sqrt(Pc3Dout.MD2eff(ndx_3D)));
                    yrngmin = 0;
                    MDlabl = MDlabel;
                else
                    yplt1_MD = MD2buf(ndx_MD);
                    yplt2_MD = MS2buf(ndx_MD);
                    yplt3_MD = Pc3Dout.MD2eff(ndx_3D);
                    yrngmin = -Inf;
                    MDlabl = MD2label;
                end
                if MD_plotting_option < 0
                    yplt = [yplt1_MD yplt2_MD yplt3_MD];
                else
                    yplt = [yplt1_MD yplt2_MD];
                end
                ymn = min(yplt);
                ymx = max(yplt);
                yrng = plot_range([ymn ymx],0.10,0.35);
                yrng(yrng < yrngmin) = yrngmin;
                plot(xplt_MD,yplt1_MD,'-c','LineWidth',lnwdthick+1);
                hold on;
                plot(xplt_MD,yplt2_MD,'-k','LineWidth',lnwdthick);
                if MD_plotting_option < 0
                    plot(xplt_3D,yplt3_MD,':m','LineWidth',lnwdthick);
                end
                hold off;
                xlim(xrng); ylim(yrng);
                % xlabel(xlabl,'FontWeight',xfwt);
                ylabel(MDlabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);

                % Plot the Ncdot rate
                subplot(4,2,7);
                xplt = minRelTimes(i)+Pc3Dout.Teph/86400;
                plot(xplt,Pc3Dout.Ncdot,'-', ...
                    'LineWidth',lnwdthick,'Color',Ncdclr);
                yrng = plot_range(Pc3Dout.Ncdot,0.10,0.25);
                yrng(yrng < 0) = 0;
                ylogplot = true;
                if ylogplot
                    set(gca,'YScale','log');
                else
                    ylim(yrng);                
                end
                xlim(xrng);
                xlabel(xlabl,'FontWeight',xfwt);
                ylabl = 'Coll.Rate (s^{-1})';
                ylabel(ylabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);
                
                % Plot the Ncdot rate relative to TCA over the conjunction
                % duration
                subplot(4,2,8);
                xrng = plot_range([Pc3Dout.TaConj Pc3Dout.TbConj],0.10);
                plot(Pc3Dout.Teph,Pc3Dout.Ncdot,'-', ...
                    'LineWidth',lnwdthick,'Color',Ncdclr);
                xlim(xrng); ylim(yrng);
                xlabl = 'Time from TCA (s)';
                xlabel(xlabl,'FontWeight',xfwt);
                ylabel(ylabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);
                
                % Plot the MD curves over the conj duration
                subplot(4,2,6);
                xplt_MD = (xplt_MD-minRelTimes(i))*86400;
                xplt_3D = (xplt_3D-minRelTimes(i))*86400;
                ndx_MD = (Pc3Dout.TaConj <= xplt_MD) & (xplt_MD <= Pc3Dout.TbConj);
                ndx_3D = (Pc3Dout.TaConj <= xplt_3D) & (xplt_3D <= Pc3Dout.TbConj);
                if MD_plotting_option < 0
                    yplt = [yplt1_MD(ndx_MD) yplt2_MD(ndx_MD) yplt3_MD(ndx_3D)];
                else
                    yplt = [yplt1_MD(ndx_MD) yplt2_MD(ndx_MD) yplt3_MD(ndx_3D)];
                end
                ymn = min(yplt);
                ymx = max(yplt);
                yrng = plot_range([ymn ymx],0.10,0.35);
                yrng(yrng < yrngmin) = yrngmin;
                plot(xplt_MD,yplt1_MD,'-c','LineWidth',lnwdthick+1);
                hold on;
                plot(xplt_MD,yplt2_MD,'-k','LineWidth',lnwdthick);
                if MD_plotting_option < 0
                    plot(xplt_3D,yplt3_MD,':m','LineWidth',lnwdthick);
                end
                hold off;
                xlim(xrng); ylim(yrng);
                % xlabel(xlabl,'FontWeight',xfwt);
                ylabel(MDlabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);
                
                % Plot the relative distance over the conj duration
                subplot(4,2,4);
                xplt = (cmbRelTimes-minRelTimes(i))*86400;
                plot(xplt,CMBEPH.RD,':k','LineWidth',lnwdthin);
                hold on;
                plot(xplt(ndx_RD),CMBEPH.RD(ndx_RD),'-k','LineWidth',lnwdthick);
                % plot(minRelTimes(i),MINEPH.RD(i),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
                hold off;
                xlim(xrng);
                set(gca, 'YLimSpec', 'Tight');
                yl = ylim;
                yrng = plot_range(yl,0.10,0.35);
                ylim(yrng);
                % xlabel(xlabl,'FontWeight',xfwt);
                ylabl = 'Separation (km)';
                ylabel(ylabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);
                
                drawnow;
                
                % Save plot;
                % orient(gcf,'Portrait');
                outputFile = [priIDStr '_' secIDStr '_' TCAtrm '.png'];
                saveas(gcf,fullfile(outputPath,outputFile));

            end

        end
        
    end
    
    % All approach event data
    AllData = struct2table(out);
    
    % Only Pc >= 0 approach event data
    cutPcLevel = 0;
    ndx1 = (AllData.Pc3D >= cutPcLevel);
    ndx2 = (AllData.Pc2D > cutPcLevel) & AllData.Pc3D < 0;
    ndxreg = ndx1 | ndx2;
    RegData = AllData(ndxreg,:);
    
    % Only Pc >= sigPcLevel approach event data
    cutPcLevel = sigPcLevel;
    ndx1 = (AllData.Pc3D >= cutPcLevel);
    ndx2 = (AllData.Pc2D > cutPcLevel) & AllData.Pc3D < 0;    
    ndxsig = ndx1 | ndx2;
    SigData = AllData(ndxsig,:);
    
    % Write event tables
    
    if writing_option >= 4
        writetable(AllData,fullfile(outputPath,[outputRoot '_all.xlsx']));
    end

    if writing_option >= 3
        writetable(RegData,fullfile(outputPath,[outputRoot '_reg.xlsx']));
    end

    if writing_option >= 2
        writetable(SigData,fullfile(outputPath,[outputRoot '_sig.xlsx']));
    end
    
    
    
    
    
    
    
    
    % Display the results
    for nResults = 0:1
        if nResults == 0
            mResults = '2D-Pc';
            PcResults = AllData.Pc2D;
        else
            mResults = '3D-Nc';
            PcResults = AllData.Pc3D;
        end
        disp([mResults ' method results:']);
        ndx = ~isnan(PcResults);
        disp([' Number of Pc values not NaN = ' num2str(sum(ndx))]);
        NcTotResults = sum(PcResults(ndx));
        disp([' NcTot = ' num2str(NcTotResults,Pcfmt)]);
        lnPsResults = sum(log(1-PcResults(ndx)));
        PcMaxResults = 1-exp(lnPsResults);
        disp([' PcMax = ' num2str(PcMaxResults,Pcfmt)]);
        if SVprocessing
            ndx = ndx & AllData.InScr > 0;
            disp(['  Number of Pc values inside SV and not NaN = ' num2str(sum(ndx))]);
            NcTotResults = sum(PcResults(ndx));
            disp(['  NcTot = ' num2str(NcTotResults,Pcfmt)]);
            lnPsResults = sum(log(1-PcResults(ndx)));
            PcMaxResults = 1-exp(lnPsResults);
            disp(['  PcMax = ' num2str(PcMaxResults,Pcfmt)]);
        end
    end
    
    
    
    
    
       
    % Process results in two passes for 2D-Pc and 3D-Pc methods
        
    idx = out.Pc3D < 0; % Events that 3D-Pc failed to converge
    if any(idx)
        Pc3D_converged = false; MethodUse = '2D-Pc';
        warning('One or more 3D-Pc calculations did not converge; not plotting results');
        nmethod1 = 1; nmethod2 = 1;
    else
        Pc3D_converged = true; MethodUse = '3D-Nc';
        nmethod1 = 1; nmethod2 = 2; % Research mode
        % nmethod1 = 2; nmethod2 = 2; % Alternate mode
    end
    
    % Process one Pc estimation method at a time
    
    for nmethod=nmethod1:nmethod2
        
        if nmethod == 1
            % 2D-Pc method analysis to estimate conjunction Pc values
            PcStr = '2D-Pc';
            MethodStr = '2D-Pc';
            PcBest = out.Pc2D;
            Tpeak = out.Tpeak2D';
            reltaua = reltaua2D;
            reltaub = reltaub2D;
            Blend = out.Blend2D;
            Long = out.Long2D;
            Var = out.NminMS2;
        else
            % 3D-Nc method analysis to estimate conjunction Pc values
            PcStr = 'Pc';
            MethodStr = '3D-Nc';
            PcBest = out.Pc3D;
            Tpeak = out.Tpeak3D';
            reltaua = reltaua3D;
            reltaub = reltaub3D;
            NPDEffect = out.NPDEffect3D;            
            Blend = out.Blend3D;
            Long = out.Long3D;
            Var = out.NmaxNcdot;
        end
        
        % Earliest and latest times to consider
        if ~any(ndxsig)
            cutPcLevel = max(PcBest)/1e3;
            ndx1 = (AllData.Pc3D >= cutPcLevel);
            ndx2 = (AllData.Pc2D > cutPcLevel) & AllData.Pc3D < 0;
            ndxplt = ndx1 | ndx2;
            if any(ndxplt)
                xbeg = min(relTA(ndxplt)); xend = max(relTB(ndxplt));
            else
                xbeg = min(relTA(ndxreg)); xend = max(relTB(ndxreg));
            end
        else
            xbeg = min(relTA(ndxsig)); xend = max(relTB(ndxsig));
        end
        
        % X-axis times bracketing overall TCA and significant-Pc
        % conjunctions
        [dca,ica] = min(MINEPH.RD);
        xbeg = min(xbeg,minRelTimes(ica));
        xend = max(xend,minRelTimes(ica));
        
        % Parameters for subplots
        Msubplot = [0.10 0.10]; % margin X & Y
        Gsubplot = [0.02 0.05]; % gutter X & Y
        Nsubplot = [2 2];
        Lsubplot = (1-Msubplot-Nsubplot.*Gsubplot)./Nsubplot;
        
        % Two plotting passes, one for Full view and one for Zoom view
        for pltpass=1:2
        
            % Plot indices
            if pltpass == 1
                % Plot all times
                ndx = true(size(cmbRelTimes));
                xpd = 0;
            else
                % Indices at times near overall TCA and significant-Pc
                % conjunctions
                ndx = (xbeg <= cmbRelTimes) & (cmbRelTimes <= xend);
                xpd = 0.1;
            end

            if plotting_option > 0

                % Set up plotting
                figure('Visible','On'); clf;
                prow = 2;

                % X plotting range
                xplt = cmbRelTimes(ndx);
                xrng = plot_range(xplt,xpd);
                % xlabl = 'Ephemeris Time (days)';
                xlabl = 'Time (days)';
                
                % Plot the conjunction Pc values, using stem plot format if there
                % are sufficiently few
                % subplot(prow,2,3);
                sp_left   = Msubplot(1);
                sp_width  = Lsubplot(1);
                sp_bottom = Msubplot(2);
                sp_height = Lsubplot(2);
                sp_pos = [sp_left sp_bottom sp_width sp_height];
                subplot('Position', sp_pos);
                
                % xxx = minRelTimes(ndxreg)+Tpeak(ndxreg)/86400;
                xxx = minRelTimes(ndxreg);
                yyy = PcBest(ndxreg);
                logmaxyyy = log10(max(yyy));
                ymx = min(1,ceil(logmaxyyy+0.1));
                ymn = min(log10(pltPcLevel),floor(logmaxyyy-3));
                yrng = 10.^plot_range([ymn ymx],0);
                % Use stems if there are few enough
                Nstem = sum(yrng(1) < yyy & yyy < yrng(2));
                if (Nstem < 10)
                    stemlsty = '-';
                    stemwid = 0.5;
                    stemsiz = msiz+2;
                elseif (Nstem < 25)
                    stemlsty = '-';
                    stemwid = 0.5;
                    stemsiz = msiz+1;
                elseif (Nstem < 100)
                    stemlsty = '-';
                    stemwid = 0.5;
                    stemsiz = msiz;
                else
                    stemlsty = 'none';
                    stemwid = 0.25;
                    stemsiz = msiz;
                end
                
                plot(NaN,NaN); hold on;
                
%                 % Plot zeros
%                 ndx = (yyy == 0);
%                 zzz = yyy(ndx)';
%                 zzz(zzz < yrng(1)) = yrng(1);
%                 clr = 'k';
%                 mrkr = 'v';
%                 stem(xxx(ndx), zzz, 'Color', clr, 'LineStyle', stemlsty, ...
%                     'LineWidth', stemwid, 'Marker', mrkr, ...
%                     'MarkerEdgeColor',clr,'MarkerFaceColor',clr, 'MarkerSize', stemsiz);
%                 hold on;
%                 % Plot tiny ones as downward triangles
%                 ndx = (yyy > 0) & (yyy <= yrng(1));
%                 zzz = yyy(ndx)';
%                 zzz(zzz < yrng(1)) = yrng(1);
%                 clr = drkgrn;
%                 mrkr = 'v';
%                 stem(xxx(ndx), zzz, 'Color', clr, 'LineStyle', stemlsty, ...
%                     'LineWidth', stemwid, 'Marker', mrkr, ...
%                     'MarkerEdgeColor',clr,'MarkerFaceColor',clr, 'MarkerSize', stemsiz);
                
                % Plot black triangles
                ndx = (yyy < yrng(1)) & (yyy < gPcLevel);
                zzz = yyy(ndx)';
                zzz(zzz < yrng(1)) = yrng(1);
                clr = 'k ';
                mrkr = 'v';
                stem(xxx(ndx), zzz, 'Color', clr, 'LineStyle', stemlsty, ...
                    'LineWidth', stemwid, 'Marker', mrkr, ...
                    'MarkerEdgeColor',clr,'MarkerFaceColor',clr, 'MarkerSize', stemsiz);
                % Plot black diamonds
                ndx = (yyy >= yrng(1)) & (yyy < gPcLevel);
                zzz = yyy(ndx)';
                clr = 'k';
                mrkr = 'd';
                stem(xxx(ndx), zzz, 'Color', clr, 'LineStyle', stemlsty, ...
                    'LineWidth', stemwid, 'Marker', mrkr, ...
                    'MarkerEdgeColor',clr,'MarkerFaceColor',clr, 'MarkerSize', stemsiz);
                % Plot green diamonds
                ndx = (yyy >= gPcLevel) & (yyy < yPcLevel);
                zzz = yyy(ndx)';
                clr = drkgrn;
                mrkr = 'd';
                stem(xxx(ndx), zzz, 'Color', clr, 'LineStyle', stemlsty, ...
                    'LineWidth', stemwid, 'Marker', mrkr, ...
                    'MarkerEdgeColor',clr,'MarkerFaceColor',clr, 'MarkerSize', stemsiz);        
                % Plot yellow diamonds
                ndx = (yyy >= yPcLevel) & (yyy < rPcLevel);
                zzz = yyy(ndx)';
                clr = gold;
                mrkr = 'd';
                stem(xxx(ndx), zzz, 'Color', clr, 'LineStyle', stemlsty, ...
                    'LineWidth', stemwid, 'Marker', mrkr, ...
                    'MarkerEdgeColor',clr,'MarkerFaceColor',clr, 'MarkerSize', stemsiz);
                % Plot red diamonds
                ndx = (yyy >= rPcLevel);
                zzz = yyy(ndx)';
                clr = 'r';
                mrkr = 'd';
                stem(xxx(ndx), zzz, 'Color', clr, 'LineStyle', stemlsty, ...
                    'LineWidth', stemwid, 'Marker', mrkr, ...
                    'MarkerEdgeColor',clr,'MarkerFaceColor',clr, 'MarkerSize', stemsiz);
                hold off;
                if SVprocessing
                    % Mark stems in screened set
                    hold on;
                    iii = round(out.InScr) ~= 0;
                    PcClip = PcBest(iii);
                    PcClip(PcClip < yrng(1)) = yrng(1);
                    plot(minRelTimes(iii),PcClip, ...
                        'LineStyle','none','LineWidth',0.5, ...
                        'Marker','.','MarkerSize',max(1,stemsiz), ...
                        'MarkerFaceColor','w','MarkerEdgeColor','w');
                    
                    hold off;
                end
                set(gca,'YScale','log');
                if prow == 2
                    xlabel(xlabl,'FontWeight',xfwt);
                end
                ylabl = 'Conjunction Pc';
                ylabel(ylabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);
                xlim(xrng); ylim(yrng);
                
%                 yt = get(gca, 'YTick');                
%                 if numel(yt) < 2
%                     ytkvct = 10.^([ymn ymx]);
%                     set(gca, 'YTick', ytkvct);        
%                 end

                if xrng(1) <= 0 && xrng(1) > -1 && ...
                   xrng(2) >= 7 && xrng(2) < 11
                    xtcks = 0:floor(xrng(2));
                    xticks(xtcks);
                else
                    xtcks = [];
                end

                latpars = [];
                latpars.NtickMin = 3;
                latpars.NtickMax = 6;
                if numel(yticklabels) < latpars.NtickMin
                    lat = LogAxisTicks(yrng,latpars);
                    idx = ~strcmpi(lat.TickLabels,'');
                    yticks(lat.Ticks(idx)); yticklabels(lat.TickLabels(idx));
                end
                grid on;
                set(gca, ...
                    'GridLineStyle','-','GridColor',gridmajcolor, ...
                    'MinorGridLineStyle','-','MinorGridColor',gridmincolor); 

            end

            % Calculate the cumulative Pc and Nc values
            PcumMax = PcBest; PcumMax(isnan(PcBest)) = 0;
            PcumMin = PcumMax; NcumEnc = PcumMax;
            for i = 2:NminTimes
                m = i-1;
                NcumEnc(i) = NcumEnc(i)+NcumEnc(m);
                PcumMin(i) = max(PcumMax(i),PcumMin(m));
                PcumMax(i) = min(1,-PcumMax(i)*PcumMax(m) + PcumMax(i) + PcumMax(m));
                % NcumEnc(i) = NcumEnc(i)+NcumEnc(i-1);
                % PcumMax(i) = 1-(1-PcumMax(i))*(1-PcumMax(i-1));
                % if (PcumMax(i) == 0) && (NcumEnc(i) > 0)
                %     PcumMax(i) = NcumEnc(i);
                % end
                % PcumMin(i) = max(PcBest(i),PcumMin(i-1));            
            end
            PcumEnd = PcumMax(end);
            if (nmethod == 1)
                % Ncum approx for Pc2D
                xxx = reltaub(ndxreg);
                yyy = NcumEnc(ndxreg)';
                Nreg = sum(ndxreg);
                if Nreg == 1
                    xbb = reltaua(ndxreg); ybb = 0;
                else
                    xbb = reltaua(ndxreg); ybb = cat(2,0,yyy(1:Nreg-1));
                end
                xxx = cat(2,xxx,xbb);
                yyy = cat(2,yyy,ybb);
                [xxx,nnn] = sort(xxx); yyy = yyy(nnn);
                if (xxx(1) > xrng(1))
                    xxx = cat(2,xrng(1),xxx);
                    yyy = cat(2,0,yyy);
                end
                if (xxx(end) < xrng(2))
                    xxx = cat(2,xxx,xrng(2));
                    yyy = cat(2,yyy,yyy(end));
                end
            else
                % Ncum for Pc3D
                xxx = Tcum;
                yyy = Ncum;
                if (Tcum(1) > xrng(1))
                    xxx = cat(2,xrng(1),xxx);
                    yyy = cat(2,0,yyy);
                end
                if (Tcum(end) < xrng(2))
                    xxx = cat(2,xxx,xrng(2));
                    yyy = cat(2,yyy,Ncum(end));
                end
            end
            NcumEnd = yyy(end);

            if plotting_option > 0

                % Plot the cumulative risk curves: Ncum, PcumMin, and PcumMax
                % subplot(prow,2,1);
                % sp_left   = Msubplot(1);
                % sp_width  = Lsubplot(1);
                sp_bottom = Msubplot(2)+Gsubplot(2)+Lsubplot(2);
                % sp_height = Lsubplot(2);
                sp_pos = [sp_left sp_bottom sp_width sp_height];
                subplot('Position', sp_pos);
                ymx = 1.05*max(yyy);
                % Plot Ncum
                plot(xxx,yyy,':','LineWidth',lnwdthin,'Color','k');
                hold on;
                % Plot PcumMin and PcumMax
                xxx = reltaub(ndxreg);
                yyy = PcumMax(ndxreg)'; 
                zzz = PcumMin(ndxreg)';
                Nreg = sum(ndxreg);
                if Nreg == 1
                    xbb = reltaua(ndxreg); ybb = 0; zbb = 0;
                else
                    xbb = reltaua(ndxreg);
                    ybb = cat(2,0,yyy(1:Nreg-1));
                    zbb = cat(2,0,zzz(1:Nreg-1));
                end
                xxx = cat(2,xxx,xbb);
                yyy = cat(2,yyy,ybb);
                zzz = cat(2,zzz,zbb);
                [xxx,nnn] = sort(xxx); yyy = yyy(nnn); zzz = zzz(nnn);
                if (xxx(1) > xrng(1))
                    xxx = cat(2,xrng(1),xxx);
                    yyy = cat(2,0,yyy);
                    zzz = cat(2,0,zzz);
                end
                if (xxx(end) < xrng(2))
                    xxx = cat(2,xxx,xrng(2));
                    yyy = cat(2,yyy,yyy(end));
                    zzz = cat(2,zzz,zzz(end));
                end
                ymx = max(ymx,1.05*max(yyy));
                yrng = [0 ymx];
                % PcumMin (for either 2D or 3D)
                plot(xxx,zzz,'--','LineWidth',lnwdthin,'Color','k');
                % PcumMax (for either 2D or 3D)            
                plot(xxx,yyy,'-','LineWidth',lnwdthin,'Color','k');
                hold off;
                % xlabel(xlabl,'FontWeight',xfwt);
                ylabl = 'Cumulative Pc';
                ylabel(ylabl,'FontWeight',yfwt);
                set(gca,'FontWeight',afwt);
                set(gca,'LineWidth',alwd);
                grid on;
                set(gca, ...
                    'GridLineStyle','-','GridColor',gridmajcolor, ...
                    'MinorGridLineStyle','-','MinorGridColor',gridmincolor); 
                xlim(xrng); ylim(yrng);
                if ~isempty(xtcks); xticks(xtcks); end
                set(gca,'Xticklabel',[]);

                % Write information to plot

                subplot(16,2,[30 32]);
                sp_left   = Msubplot(1)+Gsubplot(1)+Lsubplot(1);
                sp_bottom = Gsubplot(2);
                sp_height = Msubplot(2);
                sp_pos = [sp_left sp_bottom sp_width sp_height];
                subplot('Position', sp_pos);
                
                % plot([NaN NaN],[NaN NaN]);
                axis off;

                % Set up title structure for information
                titl = [];

                ttl = '---- Multi-conjunction Risk Analysis ----';
                titl = cat(1,titl,{ttl});
                ttl = ' ';
                titl = cat(1,titl,{ttl});

                % Write original CDM analysis creation time
                TCRstr = [TCRIDStr(1 :4 ) '-' TCRIDStr(5 :6 ) '-' TCRIDStr(7 :8 ) ' ' ... 
                          TCRIDStr(10:11) ':' TCRIDStr(12:13) ':' TCRIDStr(14:15)];
                ttl = ['Ephemeris creation time: ' TCRstr];
                titl = cat(1,titl,{ttl});

                % Write VCM file names
                if isempty(priVCMFile)
                    VCMstr = ' -- file not found --';
                else
                    VCMstr = strrep(priVCMFile,'_','\_');
                end        
                ttl = ['Primary: ' VCMstr];
                titl = cat(1,titl,{ttl});
                if isempty(secVCMFile)
                    VCMstr = ' -- file not found --';
                else
                    VCMstr = strrep(secVCMFile,'_','\_');
                end
                ttl = ['Secondary: ' VCMstr];
                titl = cat(1,titl,{ttl});

                % Write propagation time info
                ttl = ['Ephemeris start: ' dsFullBeg];
                titl = cat(1,titl,{ttl});
                ttl = ['Ephemeris end:   ' dsFullEnd];
                titl = cat(1,titl,{ttl});

                % Write the HBR
                ttl = ['Combined hard-body radius = ' num2str(HBR*1e3) ' m'];
                titl = cat(1,titl,{ttl});

                ttl = ' ';
                titl = cat(1,titl,{ttl});

                % Write the number of approaches information        
                NMDcut = sum(out.MS <= MScut);
                NPccut = sum(PcBest >= sigPcLevel);
                ttl = ['Total approach events detected: ' num2str(NminTimes)];
                titl = cat(1,titl,{ttl});
                ttl = [' (' num2str(NMDcut) ' with N{\sigma} \leq ' num2str(MScut) ...
                     ';  ' num2str(NPccut) ' with ' PcStr ' \geq ' ...
                       smart_exp_format(sigPcLevel) ')'];
                titl = cat(1,titl,{ttl});

                % Write closest approach info
                if dca < 1
                    dcastr = [num2str(dca*1e3,'%0.1f') ' m'];
                else
                    if dca < 10
                        dcastr = [num2str(dca,'%0.3f') ' km'];
                    else
                        dcastr = [num2str(dca,'%0.2f') ' km'];
                    end
                end
                ttl = ['Min conj. miss = ' dcastr ' @ ' convert_dn_to_ds(MINEPH.T(ica)+dnReference)];
                titl = cat(1,titl,{ttl});

                % Write highest Pc info
                [pmx,imx] = max(PcBest);
                Pcstr = num2str(pmx,'%0.3e');
                TPKdn = MINEPH.T(imx)+dnReference+Tpeak(imx)/86400;
                ttl = ['Max conj. Pc = ' Pcstr ' @ ' convert_dn_to_ds(TPKdn)];
                % if pmx >= rPcLevel
                %     txtclr = '{1 0 0}';
                % elseif pmx >= yPcLevel
                %     txtclr = '{1 0.667 0}';
                % else
                %     txtclr = '{0 0.667 0}';
                % end
                % ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
                titl = cat(1,titl,{ttl});

                ttl = ' ';
                titl = cat(1,titl,{ttl});

                ttl = ['Collision risks estimated using ' MethodStr ' method:'];
                titl = cat(1,titl,{ttl});

                ttl = ['Expected number of collisions: Nc = ' num2str(NcumEnd,'%0.3e')];        
                if NcumEnd >= rPcLevel
                    txtclr = '{1 0 0}';
                elseif NcumEnd >= yPcLevel
                    txtclr = '{1 0.667 0}';
                else
                    txtclr = '{0 0.667 0}';
                end
                ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
                titl = cat(1,titl,{ttl});

                PcumEnd = min([PcumEnd,NcumEnd,1]);
                if (NcumEnd > 0) && (PcumEnd == 0)
                    PcumEnd = min(1,NcumEnd);
                end
                PcBestPos = PcBest > 0; NPcBestPos = sum(PcBestPos);
                if NPcBestPos <= 1
                    opr = '=';
                elseif NPcBestPos > 1
                    opr = '\leq';
                end

                ttl = ['Accumulated collision probability: ' PcStr ...
                    ' ' opr ' ' num2str(PcumEnd,'%0.3e')];
                if PcumEnd >= rPcLevel
                    txtclr = '{1 0 0}';
                elseif PcumEnd >= yPcLevel
                    txtclr = '{1 0.667 0}';
                else
                    txtclr = '{0 0.667 0}';
                end
                ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
                titl = cat(1,titl,{ttl});

                ttl = ' ';
                titl = cat(1,titl,{ttl});

                ttl = [MethodStr ' analysis method usage violations'];
                titl = cat(1,titl,{ttl});
                ttl = ['(for events with ' PcStr ' \geq ' ...
                    smart_exp_format(sigPcLevel) '):'];
                titl = cat(1,titl,{ttl});

            else

                titl = [];

            end

            % Calculate the Pc method usage violations

            % Count violations only among significant events
            ndx0 = PcBest >= sigPcLevel; Nndx0 = sum(ndx0);

            % Initialize usage violation flag
            usage_violation = false;

            % For 3D-Pc, impose the effect of NPD TCA equinoctial covariances
            % as a potential violation
            if (nmethod == 2)
                ndx = ndx0 & (NPDEffect > NPDEffectCut); NNPDEffect = sum(ndx);
                ttl = ['Affected by NPD covariances: ' num2str(NNPDEffect) ' of ' num2str(Nndx0)];
                if (NNPDEffect > 0)
                    ttl = cat(2,ttl,'{\color{red} <violation>}');
                    usage_violation = true;
                else
                    ttl = cat(2,ttl,'{\color{blue} <no violation>}');
                end
                titl = cat(1,titl,{ttl});
            end

            % Combine blended and extended conjunction indicators into one
            % "long-duration" conjunction indicator
            ndx = ndx0 & (Blend > BlendCut | Long > LongCut); NExtend = sum(ndx);
            ttl = ['Extended or blended in time: ' num2str(NExtend) ' of ' num2str(Nndx0)];
            if (NExtend > 0)
                % Extended events are both 2D-Pc and 3D-Nc method violations
                ttl = cat(2,ttl,'{\color{red} <violation>}');
                usage_violation = true;
            else
                ttl = cat(2,ttl,'{\color{blue} <no violation>}');
            end
            titl = cat(1,titl,{ttl});

            % For 2D-Pc, impose multi-peaked events as a usage violation
            if (nmethod == 1)
                ndx = ndx0 & (Var > 1); NVar = sum(ndx);
                ttl = ['Multiple Pc-rate peaks: ' num2str(NVar) ' of ' num2str(Nndx0)];
                if (NVar > 0)
                    % Variable, multi-peaked events are 2D-Pc method violations,
                    % but not 3D-Nc method violations
                    ttl = cat(2,ttl,'{\color{red} <violation>}');
                    usage_violation = true;
                else
                    ttl = cat(2,ttl,'{\color{blue} <no violation>}');
                end
                titl = cat(1,titl,{ttl});
            end

            lPcLevel =  10^((log10(gPcLevel)+log10(yPcLevel))/2);

            auxstr = ' '; txtclr = '{0 0 1}'; BFMCpriority = 'none';
            if nmethod == 1 && Pc3D_converged
                % 2D-Pc method used but 3D-Nc method known to have converged
                recstr = ['Use ' MethodUse ' method analysis recommendation'];
                if usage_violation
                    txtclr = '{1 0 0}';
                end
            else
                if usage_violation
                    % Recommend BFMC if there is an estimation method usage violation
                    if (PcumEnd <= lPcLevel)
                        BFMCpriority = 'low';
                        txtclr = '{0 0.667 0}';
                    elseif (PcumEnd <= yPcLevel)
                        BFMCpriority = 'medium';
                        txtclr = '{1 0.667 0}';
                    else
                        BFMCpriority = 'high';
                        txtclr = '{1 0 0}';
                    end
                    recstr = ['Estimate cumulative Pc using BFMC (' BFMCpriority ' priority)'];
                else
                    recstr = ['Use ' MethodUse ' method collision probability estimates'];
                end
                if nmethod == 1
                    auxstr = '(Notify analysis team of 3D-Nc method convergence failure)';
                end            
            end

            recstr0 = recstr; auxstr0 = auxstr;

            % Complete the plotting

            if plotting_option > 0

                % Skip recommendation until reconciling with OCMDBprocess.m
                % recstr = ['{\color[rgb]' txtclr recstr '}']; %#ok<AGROW>
                % auxstr = ['{\color[rgb]' txtclr auxstr '}']; %#ok<AGROW>
                % 
                % ttl = ' ';
                % titl = cat(1,titl,{ttl});
                % ttl = '--- ANALYSIS RECOMMENDATION ---';
                % titl = cat(1,titl,{ttl});
                % 
                % ttl = recstr;
                % titl = cat(1,titl,{ttl});
                % 
                % ttl = auxstr;
                % titl = cat(1,titl,{ttl});
                
                ttl = ' ';
                titl = cat(1,titl,{ttl});
                ttl = ' ';
                titl = cat(1,titl,{ttl});
                % ttl = ' ';
                % titl = cat(1,titl,{ttl});
                % ttl = ' ';
                % titl = cat(1,titl,{ttl});
                % ttl = ' ';
                % titl = cat(1,titl,{ttl});
                % ttl = ' ';
                % titl = cat(1,titl,{ttl});

                % Write assembled title information
                htit = title(titl,'FontWeight',tfwt);
                hext = htit.Extent;
                if hext(3) < 0.8 || hext(3) > 1
                    htit.FontSize = htit.FontSize/hext(3);
                end
                
                drawnow;

                % Save plot;
                % orient(gcf,'Portrait');
                if pltpass == 1
                    pltfile = fullfile(outputPath,[outputRoot '_' ...
                       strrep(MethodStr,'-','') '_Full.png']);
                else
                    pltfile = fullfile(outputPath,[outputRoot '_' ...
                       strrep(MethodStr,'-','') '_Zoom.png']);
                end
                saveas(gcf,pltfile);

            end
            
        end
            
        % Define analysis structure
        An.Primary = {priIDStr};
        An.Secondary = {secIDStr};
        An.PcCumulative = {PcumEnd};
        An.PcMethod = {MethodStr};
        An.CATot = NminTimes;
        An.CASig = sum(PcBest >= sigPcLevel);
        if strcmpi(BFMCpriority,'none')
            An.BFMCRecommended = 'no';
        else
            An.BFMCRecommended = 'yes';
        end
        An.BFMCPriority = {BFMCpriority};
        An.BFMCStart = {[' ' dsFullBeg]};
        An.BFMCStop = {[' ' dsFullEnd]};
        if isempty(strtrim(auxstr0))
            An.Recommendation = {[recstr0 '.']};
        else
            An.Recommendation = {[recstr0 '. ' auxstr0 '.']};            
        end
        
        % Output analysis structure is final method processed
        Analysis = struct2table(An);
        
        if writing_option > 0
            writetable(Analysis,fullfile(outputPath,[outputRoot '_' ...
                strrep(MethodStr,'-','') '.xlsx']));
        end

    end
    
    return;

end

% =========================================================================

function [EphemSatID,EphemFile] = getEphemFile(SatID,inputPath)

% Get a satellite's ephemeris file name, checking for 5 or 9 digit
% satellite ID number variants.  Return empty set if no file found.

% Try the default ephemeris file name

EphemSatID = SatID;
EphemFile = fullfile(inputPath,[SatID '.eci']);
EphemFileExists = exist(EphemFile,'file');

if EphemFileExists
    return;
end

% If default file doesn't exist, then try the 5 or 9 digit name variants,
% if possible

EphemSatID = [];
EphemFile  = [];

NumID = str2double(SatID);
if isempty(NumID) || isnan(NumID) || (NumID < 0)
    % SatID does not define a usable number, so no variant file name can be
    % constructed
    return;
end

LenID = length(SatID);

if LenID == 5
    % Try 9 digit variant
    SatIDAlt = ['0000' SatID];
    EphemFileAlt = fullfile(inputPath,[SatIDAlt '.eci']);
    if exist(EphemFileAlt,'file')
        EphemSatID = SatIDAlt;
        EphemFile  = EphemFileAlt;
    end
elseif (LenID == 9) && strcmpi(SatID(1:4),'0000')
    % Try 5 digit variant
    SatIDAlt = SatID(5:9);
    EphemFileAlt = fullfile(inputPath,[SatIDAlt '.eci']);
    if exist(EphemFileAlt,'file')
        EphemSatID = SatIDAlt;
        EphemFile  = EphemFileAlt;
    end
end

return;
end

% =========================================================================

function [A] = readFortranAsTable(fileName, numElems, formatSpec, colNames)
    % Read a FORTRAN-produced data file into a Matlab table
    fileData = fileread(fileName);
    fileData = strrep(fileData,'D','E');
    sizeA = [numElems Inf];
    A = sscanf(fileData,formatSpec,sizeA);
    A = array2table(A','VariableNames',colNames);
end

% =========================================================================

function [VCMFile] = getVCMFile(SatID,inputPath)

% Get a satellite's VCM file name, checking for 5 or 9 digit
% satellite ID number variants.  Return empty set if no file found.

% Initialize output

VCMFile  = [];

% Search using input SatID

template = fullfile(inputPath,[SatID '*.vcm']);
dl = dir(template);
if ~isempty(dl)
    if numel(dl) > 1
        warning('More than one candidate VCM file found; using first');
        VCMFile = [dl(1).name ' (among ' numel(dl) ' candidates ']; 
    else
        VCMFile = dl.name;
    end
    return;
end

% If default file template doesn't work, then try the 5 or 9 digit name
% variants, if possible.

NumID = str2double(SatID);
if isempty(NumID) || isnan(NumID) || (NumID < 0)
    % SatID does not define a usable number, so no variant file name can be
    % constructed
    return;
end

LenID = length(SatID);

if LenID == 5
    % Try 9 digit variant
    SatIDAlt = ['0000' SatID];
    template = fullfile(inputPath,[SatIDAlt '*.vcm']);
    dl = dir(template);
elseif (LenID == 9) && strcmpi(SatID(1:4),'0000')
    % Try 5 digit variant
    SatIDAlt = SatID(5:9);
    template = fullfile(inputPath,[SatIDAlt '*.vcm']);
    dl = dir(template);
end

if ~isempty(dl)
    if numel(dl) > 1
        warning('More than one candidate VCM file found; using first');
        VCMFile = [dl(1).name ' (among ' numel(dl) ' candidates ']; 
    else
        VCMFile = dl.name;
    end
end

return;
end

% =========================================================================

function [S,M2,M2HBR] = calc_sep_maha_distances(EPH,HBR,Lclip)
% Calculate separation and Mahalanobis distances from ephemeris
N = numel(EPH.T);
r = EPH.X2(1:3,:)-EPH.X1(1:3,:);
S = sqrt(sum(r.*r,1));
if nargout < 2
    return
end
M2 = NaN(size(S)); M2HBR = M2;
for n=1:N
    A = EPH.P1(1:3,1:3,n)+EPH.P2(1:3,1:3,n);
    [M2(n),M2HBR(n)] = calc_maha_distances(r(:,n),A,HBR,Lclip);
end
return
end

% =========================================================================

function [M2,M2HBR] = calc_maha_distances(r,A,HBR,Lclip)
% Calculate Mahalanobis distances from rel.pos. vector and rel.pos
% covariance matrix.
[~,~,~,~,~,~,Ainv] = CovRemEigValClip(A(1:3,1:3),Lclip);
Air = Ainv * r;
% MD2 for center of collision sphere w/ radius HBR; this should be
% nonnegative
M2 = r' * Air;
% Min. MD2 over collision sphere, to 1st order accuracy in HBR/|r|; this
% approximation can yield negative estimates
M2HBR = M2 - 2*HBR*norm(Air);
return
end

% =========================================================================

function MS2 = interp_MS2(T,EPH,HBR,Lclip,verbose)

% Calculate MS2 at times T by interpolating an ephemeris

N = numel(EPH.TREL);
if any(T < EPH.TREL(1)) || any(T > EPH.TREL(N))
    error('Time(s) out of ephemeris bounds');
elseif (N < 5)
    error('Too few ephemeris points');
end

M = numel(T);
MS2 = NaN(size(T));

for m=1:M
        
    % Find best five bracketing eph points for interpolation
    [~,n] = min(abs(T(m)-EPH.TREL));
    n1 = n-2;
    n2 = n+2;
    if (n1 < 1)
        n1 = 1; n2 = 5;
    elseif (n2 > N)
        n1 = N-4; n2 = N;
    end

    % Interpolate state and covariance
    i = n1:n2;
    [r1,~,P1] = StateCovInterp(T(m), ...
        EPH.TREL(i), EPH.X1(1:3,i)', EPH.X1(4:6,i)', EPH.P1(:,:,i));
    [r2,~,P2] = StateCovInterp(T(m), ...
        EPH.TREL(i), EPH.X2(1:3,i)', EPH.X2(4:6,i)', EPH.P2(:,:,i));
    
    % If appropriate, interpolate states with the alternate five bracketing
    % points, and create weighted average states
    if (T(m) ~= EPH.TREL(n))
        % If not at eph point, generate the alternate five bracketing points
        if T(m) > EPH.TREL(n)
            k = n+1;
        else
            k = n-1;
        end
        k1 = k-2;
        k2 = k+2;
        if (k1 < 1)
            k1 = 1; k2 = 5;
        elseif (k2 > N)
            k1 = N-4; k2 = N;
        end
        if (k1 ~= n1)
            % Alternate bracketing points are different than the original,
            % so perform the interpolation
            i = k1:k2;
            [a1,~,Q1] = StateCovInterp(T(m), ...
                EPH.TREL(i), EPH.X1(1:3,i)', EPH.X1(4:6,i)', EPH.P1(:,:,i));
            [a2,~,Q2] = StateCovInterp(T(m), ...
                EPH.TREL(i), EPH.X2(1:3,i)', EPH.X2(4:6,i)', EPH.P2(:,:,i));
            % Weight the two alternate sets of states by time differences
            w = (T(m)-EPH.TREL(n)) / (EPH.TREL(k)-EPH.TREL(n));
            if (w <= 0) || (w >= 1)
                error('Invalid weight');
            end
            omw = 1-w;
            r1 = omw*r1+w*a1; P1 = omw*P1+w*Q1;
            r2 = omw*r2+w*a2; P2 = omw*P2+w*Q2;
        end
    end
    
    % Calculate the MS2 value
    r = (r2-r1)';
    A = P1(1:3,1:3)+P2(1:3,1:3);
    [~,MS2(m)] = calc_maha_distances(r,A,HBR,Lclip);
    
    if verbose
        disp(num2str(n1:n2,'%2i  '));
        disp(num2str(r1,'%0.16e  '));
        disp(num2str(P1(1:3,1:3),'%0.16e  '));
        disp(num2str(r2,'%0.16e  '));
        disp(num2str(P2(1:3,1:3),'%0.16e  '));
        disp(num2str(r','%0.16e  '));
        disp(num2str(A,'%0.16e  '));
    end
    
end
return
end

% =========================================================================

function [Tcum] = calc_cum_times(T)
N = numel(T);
M = N-1;
Tcum = [0.5*(T(1:M)+T(2:N)) T(N)];
% Tcum = [0.5*(T(1:M)+T(2:N)) T(N)+0.5*(T(N)-T(M))];
return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
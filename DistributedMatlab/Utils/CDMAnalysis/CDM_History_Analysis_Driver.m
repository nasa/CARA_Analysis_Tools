%% CDM File Analysis Driver
% 
% Purpose: This code is intended to act as a stand alone Driver to parse and process
% an input CDM txt file and report Probability of Collision estimates using
% multiple Methods

% Clear Existing variables from workspace
clear

%% Set Up File Paths
p = mfilename('fullpath');
[filepath,~,~] = fileparts(p);
cd([filepath filesep '..' filesep '..']);
addpath(genpath(pwd))
ParentPath = pwd;

%% Constants
FiguresFolder   = '';
% FiguresFolder   = 'D:\TemporaryFigures\';
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
         
         
RICCovScaling = eye(6); RICCovScaling(3,3) =3;
RICCovScaling = eye(6);

% Set up Color Schemes
RGBTriplets = [0 0.4470 0.7410
               0.8500 0.3250 0.0980
               0.9290 0.6940 0.1250
               0.4940 0.1840 0.5560
               0.4660 0.6740 0.1880
               0.3010 0.7450 0.9330
               0.6350 0.0780 0.1840];


%% Allow User to Select Multiple CDMs and Set Inputs
% Retrieve File
[FileName,PathName] = uigetfile({'*.txt;*.cdm';'*.*'},'Please Select CDM Input File(s)','MultiSelect','on');
if ~iscell(FileName) && FileName==0
    fprintf('\nNo file selected, program will now exit\n')
    return
end

% Parse CDM Input Files
CDMfile = cell(length(FileName),1);
for i=1:length(FileName)
    CDMfile{i} = fullfile(PathName, FileName{i});
end
% [cdmhead, cdmobj] = read_cdm(CDMfile);

% Get Hard Body Radius (m) if required
HBR         = 20;
answer      = inputdlg('Input Hard Body Radius (HBR) in Meters','HBR',1,{num2str(HBR)});
if isempty(answer)
    fprintf('\nNo Hard Body Radius Input, program will now exit\n')
    return
end
HBR         = str2double(answer{1});

% Set Output folder for generated Figures if required
if isempty(FiguresFolder)
    FiguresFolder = uigetdir('','Please Select A Directory to save figures to (optional)');
    if FiguresFolder==0
        fprintf('\nNo output folder selected, program will continue though plots will not be saved\n')
    end
end

% get filename of obs file if available
[ObsFile, ObsPath] = uigetfile('*.*','Please Select A Relevant Observation record file if available (optional)');
if ObsFile==0
    fprintf('\nNo Observation File selected, program will not plot actual tracks after sensor coverage tool\n')
else
    ObsFile = fullfile(ObsPath,ObsFile);
end

% get filenames of solar_flux data files
[FileName,PathName] = uigetfile('*.*','Please Select solar_flux Input File(s)','MultiSelect','on');
if ~iscell(FileName) && FileName==0
    fprintf('\nNo Solar Flux File Selected, program will not produce predictive plots\n')
    SolarFluxfile = {};
else
    SolarFluxfile = cell(length(FileName),1);
    for i=1:length(FileName)
        SolarFluxfile{i} = fullfile(PathName, FileName{i});
    end
end
%% Get Space Weather Data
% retrieve from online
SpaceWeatherFile = fullfile(filepath,'SW.txt');
SpaceWeatherFile = websave(SpaceWeatherFile,'https://celestrak.com/SpaceData/SW-Last5Years.txt');

% parse data
SWData = parseSWData(SpaceWeatherFile);

%% Execute Pc Omnibus Tool
MessageTimes = zeros(length(CDMfile),1);
for i = 1:length(CDMfile)
    cdmName = strsplit(CDMfile{i},{'.',filesep});
    cdmName = cdmName(end-1);
    [PcStructTemp]   = Pc_Omnibus(CDMfile{i},HBR,fullfile(FiguresFolder,cdmName{1}),[],false);
    PcStructTemp.CDM.OBJECT2_C_RTN = RICCovScaling*PcStructTemp.CDM.OBJECT2_C_RTN*RICCovScaling;
    if i==1
        PcStructOut = PcStructTemp;
    else
        PcStructOut = [PcStructOut; PcStructTemp];
    end
    close all
    MessageTimes(i) = datenum(PcStructTemp.CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS');
end

% reorder by message time
[~,idx]         = sort(MessageTimes);
MessageTimes    = MessageTimes(idx);
PcStructOut     = PcStructOut(idx);

if isfield(PcStructOut(1).CDM,'OBJECT1_OBJECT_NAME') && ~isempty(PcStructOut(1).CDM.OBJECT1_OBJECT_NAME)
    PrimaryText = [PcStructOut(1).CDM.OBJECT1_OBJECT_NAME '(' PcStructOut(1).CDM.OBJECT1_OBJECT_DESIGNATOR ')'];
else
    PrimaryText = PcStructOut(1).CDM.OBJECT1_OBJECT_DESIGNATOR;
end

%% Generate OD history plot

% Calculate Max and Minimum LUPIs
for i = 1:length(PcStructOut)
    [EDRBin] = EDRBinDetermination(PcStructOut(i).CDM.OBJECT2_SEDR);
    period = orbit_period([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z]*1000,[PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT]*1000,GM);
    period = period/60; % Convert to minutes
    [KEP] = Cart2Kep([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT],'True','Deg');
    ecc = KEP(2);
    [PcStructOut(i).maxLUPI] = MaximumLUPIDetermination(EDRBin,ecc,period);
    [PcStructOut(i).minLUPI] = MinimumLUPIDetermination(EDRBin,ecc);
end


% Load Space-Track Credentials

% Prompt user if none available

% Save credentials locally

% retrieve Current TLE history

% Generate Plot
FigurePos       = [50 50 600 600];
ColorMapDef     = [1 0 1 0.5%Magenta
                   0 0 1 0.5]; %Blue
y_limits = [0 0];
x_limits = [0 ceil(datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')+0.5)];
g = figure('Name','ODHistory',...
           'Position', FigurePos);
hold on
% plot OD spans
for i = 1:length(PcStructOut)
    epoch = datenum(strrep(PcStructOut(i).CDM.OBJECT2_TIME_LASTOB_END,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF');
    [~,CompositeScore,~]=CDM_ODQualityAssessment(PcStructOut(i).CDM,[],[]);
    if PcStructOut(i).maxLUPI > PcStructOut(i).CDM.OBJECT2_ACTUAL_OD_SPAN
        plot([epoch-PcStructOut(i).maxLUPI, epoch-PcStructOut(i).CDM.OBJECT2_ACTUAL_OD_SPAN],[epoch, epoch],'LineWidth',15,'Color',[0 1 0 0.5]); % Green Bar
        plot([epoch-PcStructOut(i).CDM.OBJECT2_ACTUAL_OD_SPAN, epoch],[epoch, epoch],'LineWidth',15,'Color',(ColorMapDef(2,:)-ColorMapDef(1,:))*CompositeScore+ColorMapDef(1,:)); % blue Bar'
    else
        plot([epoch-PcStructOut(i).CDM.OBJECT2_ACTUAL_OD_SPAN, epoch-PcStructOut(i).maxLUPI],[epoch, epoch],'LineWidth',15,'Color',[1 0 0 0.5]); % red Bar'
        plot([epoch-PcStructOut(i).maxLUPI, epoch],[epoch, epoch],'LineWidth',15,'Color',(ColorMapDef(2,:)-ColorMapDef(1,:))*CompositeScore+ColorMapDef(1,:)); % blue Bar
    end
    if i==1
        y_limits(1) = ceil(epoch-3);
        x_limits(1) = ceil(epoch-PcStructOut(i).maxLUPI-3);
    elseif i==length(PcStructOut)
        y_limits(2) = floor(epoch+3);
    end
end
% Plot Observations if available
if ObsFile ~=0
    try % Try parsing observation record file
        [ObservationStruct,TrackStruct] = ParseObservationRecordFile(ObsFile);
    catch % else try reading file as csv format
        [ObservationStruct,TrackStruct] = ParseObservationRecordFileCSV(ObsFile);
    end
    epoch = datenum(strrep(PcStructOut(end).CDM.OBJECT2_TIME_LASTOB_END,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF');
    ObsIndexes = find([TrackStruct.LastObsTime]+30/1440 >= (datenum(strrep(PcStructOut(1).CDM.OBJECT2_TIME_LASTOB_END,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF') - PcStructOut(1).CDM.OBJECT2_ACTUAL_OD_SPAN) & ...
                      [TrackStruct.FirstObsTime]-30/1440 <= epoch);
    SensorNumbers       = unique([TrackStruct(ObsIndexes).SensorNumber]);
    TrackingMarkers = {'o' ,  's',  'd',  '^',  'v'};
    MarkerUnicodes  = [9679, 9632, 9830, 9650, 9660];
    for i = 1:length(PcStructOut)
        epoch = datenum(strrep(PcStructOut(i).CDM.OBJECT2_TIME_LASTOB_END,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF');
        ObsIndexes = find([TrackStruct.LastObsTime]+30/1440 >= (epoch - PcStructOut(i).CDM.OBJECT2_ACTUAL_OD_SPAN) & ...
                          [TrackStruct.FirstObsTime]-30/1440 <= epoch);

        for j=1:length(SensorNumbers)
            PlotIndexes     = find([TrackStruct(ObsIndexes).SensorNumber] == SensorNumbers(j));
            if ~isempty(PlotIndexes)
                scatter([TrackStruct(ObsIndexes(PlotIndexes)).FirstObsTime],...
                        epoch * ones(size([TrackStruct(ObsIndexes(PlotIndexes)).FirstObsTime])),...
                        'filled',TrackingMarkers{mod(j-1,length(TrackingMarkers))+1},...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor','k')
            end
        end 
    end
    if length(SensorNumbers)<=length(TrackingMarkers)
        AxisPosition = get(gca,'Position');
        SensorLegend = {'Sensor IDs'};
        for j=1:length(SensorNumbers)
            SensorLegend = [SensorLegend;
                            {[char(MarkerUnicodes(j)) ': ' num2str(SensorNumbers(j),'%d')]}];
        end
        annotation('textbox',[AxisPosition(1) AxisPosition(2) 0.145*FigurePos(3)/600 length(SensorNumbers)*25/FigurePos(4)],...
           'String',SensorLegend,...
           'LineStyle','-',...
           'BackgroundColor','w',...
           'FontUnits','points',...
           'FontSize',10,...
           'FontName','Consolas',...
           'FaceAlpha',1);
    end
end
% plot TCA
plot([datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF') datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')],y_limits,'--r');
hold off
title({['OD History for Object: ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR]; ['\color{red}TCA: ' PcStructOut(1).CDM.TCA ' UTC']})
ylabel('OD Epoch')
datetick('y',6)
xlabel('Observation Times')
datetick('x',6)
ax=gca;
set(ax,'Ydir','reverse')
grid on
box on

xlim(x_limits);
ylim(y_limits);

% Save figure
saveas(gca,fullfile(FiguresFolder,'ODHistory.fig'));
saveas(gca,fullfile(FiguresFolder,'ODHistory.png'));

%% Generate Sensor Coverage Map/History
try
    SensorCoverageGUI = guidata(SensorCoverage);
    
    % Set STF FilePath
    set(SensorCoverageGUI.STFPath,'String',fullfile(ParentPath,'Main','SensorCoverageExternal','stfData.mat'));
    
    % Set Folder of CDM data
    set(SensorCoverageGUI.FilePathTextBox,'String',fileparts(CDMfile{1}));
    
    % Set Primary Object
    set(SensorCoverageGUI.PrimaryTextBox,'String',PcStructOut(1).CDM.OBJECT1_OBJECT_DESIGNATOR);
    
    % Set Secondary Object
    set(SensorCoverageGUI.SecondaryTextBox,'String',PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR);
    
    % Survey Data
    SensorCoverage('SurveyData_Callback',SensorCoverageGUI.SurveyData, [], SensorCoverageGUI);
    
    % Deselect Historical Graphs
    set(SensorCoverageGUI.HistGraph,'Value',0);
    
    % Find and Select First Available CDM in series
    DataTable = get(SensorCoverageGUI.GUIDataTable,'Data');
    [~,idx] = sort(datenum(DataTable(:,4),'yyyy-mm-dd HH:MM:SS'),'ascend');
    DataTable{idx(1),7} = true;
    set(SensorCoverageGUI.GUIDataTable,'Data',DataTable);
    
    % Generate Coverage Map
    SensorCoverage('SubmitOCMButton_Callback',SensorCoverageGUI.SubmitOCMButton, [], SensorCoverageGUI);
    
    % Close GUI Data
    close(SensorCoverageGUI.figure1);
    
    % Get global variable from Sensor Coverage Tool
    global GLOBALDATA
   
    % Get Sensor Coverage figure
    SensorCovFig  = get(gca);
    % Set Sensor Coverage span limits to be similar to OD History Figure,
    % but close to message times
    set(gca,'xlim',y_limits+[2 -2]);
    % Get current sensor coverage y limits
    SensorYLimits = get(gca,'ylim');
    % Add indicator lines for each CDM OD update
    hold on
    for ii = 1:length(MessageTimes)
        h(ii) = plot([MessageTimes(ii) MessageTimes(ii)],SensorYLimits,'--b','LineWidth',2);
        h(ii).DisplayName = [datestr(datenum(PcStructOut(ii).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),'yyyy-mm-dd HH:MM') ...
                        ' | ' num2str(PcStructOut(ii).PC2DBase,'%0.2E')];
    end
    hold off
    
    
    % Plot Actual Observations if available
    if ObsFile ~=0
        y_ticks = get(gca,'YTick');
        y_tick_labels = get(gca,'YTickLabel');
        try % Try parsing observation record file
            [ObservationStruct,TrackStruct] = ParseObservationRecordFile(ObsFile);
        catch % else try reading file as csv format
            [ObservationStruct,TrackStruct] = ParseObservationRecordFileCSV(ObsFile);
        end
            
        tempTrackStruct = TrackStruct([TrackStruct(:).FirstObsTime] >= y_limits(1) & [TrackStruct(:).FirstObsTime] <= y_limits(2));
        hold on
        for ii = 1:length(tempTrackStruct)
            % Get Sensor Name
            idx1 = find(GLOBALDATA.passes_by_site.SenID == tempTrackStruct(ii).SensorNumber,1,'first');
            if ~isempty(idx1)
                idx2 = find(~cellfun(@isempty,regexpi(y_tick_labels,GLOBALDATA.passes_by_site.SensorName(idx1))));
                if ~isempty(idx2)
                    tempTrackStruct(ii).SensorLabel = y_tick_labels(idx2);
                    tempTrackStruct(ii).ObsTimeStr = {datestr(tempTrackStruct(ii).FirstObsTime,'yyyy-mm-dd HH:MM:SS')};
                    scatter(tempTrackStruct(ii).FirstObsTime,y_ticks(idx2),40,'x','r','LineWidth',2)
                end
            end
        end
        hold off
        
        % Append Title to Indicate Actual Tracks
        CurrentTitle = get(gca,'title');
        CurrentTitle = CurrentTitle.String;
        % Ensure not constant Appending
        if isempty(regexpi(CurrentTitle{2},'\color{red}'))
            CurrentTitle{2} = [strtrim(CurrentTitle{2}) ' - \color{red}X \color{black}- Indicates Actual Track'];
            title(gca,CurrentTitle);
        end
    end
    Leg = legend(h,'LineWidth',1,'NumColumns',1);
    title(Leg,'       Message Time  |     Pc   ');
    
    % Save figure
    saveas(gca,fullfile(FiguresFolder,'SensorCoverage.fig'));
    saveas(gca,fullfile(FiguresFolder,'SensorCoverage.png'));
catch
    % Close GUI Data
    close(SensorCoverageGUI.figure1);
    fprintf(['\nNo tracking data available for object: ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ', Sensor coverage cannot be run for this object. Program will skip sensor coverage and now continue\n'])
end
%% Generate Individual Object RIC Uncertainty History
figure('NumberTitle','off','Name','RIC Miss/Uncertainty and Mahalanobis Update Plots','Position',[500,150,650,800])
% Mahalonobis Distance
subplot(4,1,1)
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Individual Object Mahalanobis Distance and RIC Miss Components Update Cycles'});
hold on
i=1;
SecondaryDiff = zeros(length(CDMfile),3);
SecondaryDiff(i,:) = zeros(3,1);
PrimaryDiff = zeros(length(CDMfile),3);
RICCombDiff = zeros(length(CDMfile),3);
RICUpdateDiff = zeros(length(CDMfile),6);
PrimaryDiff(i,:) = zeros(3,1);
AngleDiff = zeros(length(CDMfile),1);
MDCorr = zeros(length(CDMfile)-1,5);
clear h
for i = 2:length(CDMfile)
% This commented section was to check if variations had to do with ITRF vs
% ECI comparisons
    [r_s2,v_s2] = PosVelConvert([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z],...
                                [PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT],...
                                strrep(PcStructOut(i).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    [r_s1,v_s1] = PosVelConvert([PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z],...
                                [PcStructOut(i-1).CDM.OBJECT2_X_DOT PcStructOut(i-1).CDM.OBJECT2_Y_DOT PcStructOut(i-1).CDM.OBJECT2_Z_DOT],...
                                strrep(PcStructOut(i-1).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    SecondaryDiff(i,:) = (r_s2 ...
                          +(datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
                          *v_s2 ...
                          )-r_s1;
    RICUpdateDiff(i,1:3) = ECI2RIC(SecondaryDiff(i,:),r_s1,v_s1)*1000;
    [r_s2,v_s2] = PosVelConvert([PcStructOut(i).CDM.OBJECT1_X PcStructOut(i).CDM.OBJECT1_Y PcStructOut(i).CDM.OBJECT1_Z],...
                                [PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT],...
                                strrep(PcStructOut(i).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    [r_s1,v_s1] = PosVelConvert([PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z],...
                                [PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT],...
                                strrep(PcStructOut(i-1).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    PrimaryDiff(i,:)   = (r_s2 ...
                          +(datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
                          *v_s2 ...
                          )-r_s1;
    RICUpdateDiff(i,4:6) = ECI2RIC(PrimaryDiff(i,:),r_s1,v_s1)*1000;
%     SecondaryDiff(i,:) = ([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z] ...
%                           +(datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
%                           *[PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT] ...
%                           )-[PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z];
%     PrimaryDiff(i,:)   = ([PcStructOut(i).CDM.OBJECT1_X PcStructOut(i).CDM.OBJECT1_Y PcStructOut(i).CDM.OBJECT1_Z]...
%                           +(datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
%                           *[PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT] ...
%                           )-[PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z];
    AngleDiff(i)       =  acosd(dot([PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT],[PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT])...
                                /norm([PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT])...
                                /norm([PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT]));
%     RICUpdateDiff(i,1:3) = ECI2RIC(SecondaryDiff(i,:),[PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z],[PcStructOut(i-1).CDM.OBJECT2_X_DOT PcStructOut(i-1).CDM.OBJECT2_Y_DOT PcStructOut(i-1).CDM.OBJECT2_Z_DOT])*1000;
%     RICUpdateDiff(i,4:6) = ECI2RIC(PrimaryDiff(i,:),[PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z],[PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT])*1000;
    MDCorr(i-1,1)  = sqrt(SecondaryDiff(i,:)*1000*inv(PcStructOut(i-1).CDM.OBJECT2_C_J2K(1:3,1:3))*SecondaryDiff(i,:)'*1000);
    MDCorr(i-1,1)  = sqrt(RICUpdateDiff(i,1:3)*inv(PcStructOut(i-1).CDM.OBJECT2_C_RTN(1:3,1:3))*RICUpdateDiff(i,1:3)');
    MDCorr(i-1,2)  = sqrt(PrimaryDiff(i,:)*1000*inv(PcStructOut(i-1).CDM.OBJECT1_C_J2K(1:3,1:3))*PrimaryDiff(i,:)'*1000);
    MDCorr(i-1,2)  = sqrt(RICUpdateDiff(i,4:6)*inv(PcStructOut(i-1).CDM.OBJECT1_C_RTN(1:3,1:3))*RICUpdateDiff(i,4:6)');
    RICCombDiff(i,:)    = [PcStructOut(i).CDM.RELATIVE_POSITION_R PcStructOut(i).CDM.RELATIVE_POSITION_T PcStructOut(i).CDM.RELATIVE_POSITION_N]...
                           - [PcStructOut(i-1).CDM.RELATIVE_POSITION_R PcStructOut(i-1).CDM.RELATIVE_POSITION_T PcStructOut(i-1).CDM.RELATIVE_POSITION_N];
    MDCorr(i-1,3)  = sqrt(RICCombDiff(i,:)*inv(PcStructOut(i-1).CDM.COMBINED_C_RTN(1:3,1:3))*RICCombDiff(i,:)');
    % Velocity Mahalanobis Distances
    MDCorr(i-1,4)  = sqrt(([PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT]*1000 ...
                           -[PcStructOut(i-1).CDM.OBJECT2_X_DOT PcStructOut(i-1).CDM.OBJECT2_Y_DOT PcStructOut(i-1).CDM.OBJECT2_Z_DOT]*1000)...
                           *inv(PcStructOut(i-1).CDM.OBJECT2_C_J2K(4:6,4:6))...
                           *([PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT]*1000 ...
                           -[PcStructOut(i-1).CDM.OBJECT2_X_DOT PcStructOut(i-1).CDM.OBJECT2_Y_DOT PcStructOut(i-1).CDM.OBJECT2_Z_DOT]*1000)');
    MDCorr(i-1,5)  = sqrt(([PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT]*1000 ...
                           -[PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT]*1000)...
                           *inv(PcStructOut(i-1).CDM.OBJECT1_C_J2K(4:6,4:6))...
                           *([PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT]*1000 ...
                           -[PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT]*1000)');
%     RICCombDiff(i,:)    = [PcStructOut(i).CDM.RELATIVE_POSITION_R PcStructOut(i).CDM.RELATIVE_POSITION_T PcStructOut(i).CDM.RELATIVE_POSITION_N];
%     MDCorr(i-1,3)  = sqrt(RICCombDiff(i,:)*inv(PcStructOut(i).CDM.COMBINED_C_RTN(1:3,1:3))*RICCombDiff(i,:)');
%     MDCorr(i-1,3)  = sqrt(RICCombDiff(i,:)*inv(diag(diag(PcStructOut(i).CDM.COMBINED_C_RTN(1:3,1:3))))*RICCombDiff(i,:)');
    if i==2
        h(1) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,1),80,'rs','filled');  
        h(2) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,2),80,'bo','filled');  
%         h(3) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,3),80,'g^','filled');  
        h(1).DisplayName = 'Secondary';
        h(2).DisplayName = 'Primary';
%         h(3).DisplayName = 'Combined';
    else
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,1),80,'rs','filled');  
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,2),80,'bo','filled');  
%         scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,3),80,'g^','filled');  
    end
end
hold off
box on
grid on
ylabel({'Mahalanobis'; 'Correction'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
if max(max(MDCorr(:,1:2)))>10
    set(gca, 'YScale', 'log')
    ylim([10^floor(log10(min(min(MDCorr(:,1:2))))) 10^ceil(log10(max(max(MDCorr(:,1:2)))))])
end

legend(h);

% Radial Miss and Combined Uncertainty
subplot(4,1,2)
hold on
for i = 1:length(CDMfile)
%     errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.RELATIVE_POSITION_R,sqrt(PcStructOut(i).CDM.COMBINED_C_RTN(1,1)),'s','Color',[0 0.4470 0.7410],'MarkerSize',5,'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',1);
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,1),sqrt(PcStructOut(i).CDM.OBJECT2_C_RTN(1,1)),'s','Color','r','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,4),sqrt(PcStructOut(i).CDM.OBJECT1_C_RTN(1,1)),'o','Color','b','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
    
end
plot([floor(min(MessageTimes)) ceil(max(MessageTimes))],[0 0],'-k','LineWidth',0.5)
hold off
ylabel({'Radial Update'; 'Sequence (m)'})
box on
grid on
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')

% Intrack Miss and Combined Uncertainty
subplot(4,1,3)
hold on
for i = 1:length(CDMfile)
%     errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.RELATIVE_POSITION_T,sqrt(PcStructOut(i).CDM.COMBINED_C_RTN(2,2)),'s','Color',[0.8500 0.3250 0.0980],'MarkerSize',5,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',1);
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,2),sqrt(PcStructOut(i).CDM.OBJECT2_C_RTN(2,2)),'s','Color','r','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,5),sqrt(PcStructOut(i).CDM.OBJECT1_C_RTN(2,2)),'o','Color','b','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
end
plot([floor(min(MessageTimes)) ceil(max(MessageTimes))],[0 0],'-k','LineWidth',0.5)
hold off
box on
grid on
ylabel({'Intrack Update'; 'Sequence (m)'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')

% Crosstrack miss and combined uncertainty
subplot(4,1,4)
hold on
for i = 1:length(CDMfile)
%     errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.RELATIVE_POSITION_N,sqrt(PcStructOut(i).CDM.COMBINED_C_RTN(3,3)),'s','Color',[0.4660 0.6740 0.1880],'MarkerSize',5,'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerFaceColor',[0.4660 0.6740 0.1880],'LineWidth',1);
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,3),sqrt(PcStructOut(i).CDM.OBJECT2_C_RTN(3,3)),'s','Color','r','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r','LineWidth',1);
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,6),sqrt(PcStructOut(i).CDM.OBJECT1_C_RTN(3,3)),'o','Color','b','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b','LineWidth',1);
end
plot([floor(min(MessageTimes)) ceil(max(MessageTimes))],[0 0],'-k','LineWidth',0.5)
hold off
box on
grid on
ylabel({'CrossTrack Update'; 'Sequence (m)'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
xlabel('Date of Message Update');

saveas(gca,fullfile(FiguresFolder,'IndividualRICUncertaintyHistory.fig'));
saveas(gca,fullfile(FiguresFolder,'IndividualRICUncertaintyHistory.png'));
%% Generate Secondary Object Radial Plane History
figure('NumberTitle','off','Name','Secondary Radial Plane Update Plot','Position',[500,150,700,700])
% Mahalonobis Distance
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Secondary Object Radial Plane Update Cycles'});
hold on
clear h
% for i = 2:length(CDMfile)
%     SecondaryDiff(i,:) = ([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z] ...
%                           +(datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
%                           *[PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT] ...
%                           )-[PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z];
%     PrimaryDiff(i,:)   = ([PcStructOut(i).CDM.OBJECT1_X PcStructOut(i).CDM.OBJECT1_Y PcStructOut(i).CDM.OBJECT1_Z]...
%                           +(datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
%                           *[PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT] ...
%                           )-[PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z];
%     RICUpdateDiff(i,1:3) = ECI2RIC(SecondaryDiff(i,:),[PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z],[PcStructOut(i-1).CDM.OBJECT2_X_DOT PcStructOut(i-1).CDM.OBJECT2_Y_DOT PcStructOut(i-1).CDM.OBJECT2_Z_DOT])*1000;
%     RICUpdateDiff(i,4:6) = ECI2RIC(PrimaryDiff(i,:),[PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z],[PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT])*1000;
% end
for i=1:length(CDMfile)
    % Build 2X2 covariance matrix
    A = [PcStructOut(i).CDM.OBJECT2_CT_T PcStructOut(i).CDM.OBJECT2_CN_T;
         PcStructOut(i).CDM.OBJECT2_CN_T PcStructOut(i).CDM.OBJECT2_CN_N];
    % Get Eigen values
    V=eig(A);
%     V=sort(V,'descend');
    
%     a = sqrt(PcStructOut(i).CDM.OBJECT2_CT_T);
%     b = sqrt(PcStructOut(i).CDM.OBJECT2_CN_N);
%     ClockAngle = atan(PcStructOut(i).CDM.OBJECT2_CN_T/sqrt(PcStructOut(i).CDM.OBJECT2_CT_T)/sqrt(PcStructOut(i).CDM.OBJECT2_CN_N));
    a=V(1);
    b=V(2);
    ClockAngle = atan2(a,b);
    a=sqrt(a);
    b=sqrt(b);
    % Generate basis ellipse (parameterization)
    u  = (0:0.1:2*pi+0.1)';
    x0 = a * cos(u);
    y0 = b * sin(u);

    % Rotate ellipse by clock angle
    x = RICUpdateDiff(i,3)+x0.*cosd(ClockAngle) - y0.*sind(ClockAngle); 
    y = RICUpdateDiff(i,2)+x0.*sind(ClockAngle) + y0.*cosd(ClockAngle);
    h(i) = plot(x,y);

%     h(i) = errorbar(RICUpdateDiff(i,3),RICUpdateDiff(i,2),sqrt(PcStructOut(i).CDM.OBJECT2_CT_T),sqrt(PcStructOut(i).CDM.OBJECT2_CT_T),sqrt(PcStructOut(i).CDM.OBJECT2_CN_N),sqrt(PcStructOut(i).CDM.OBJECT2_CN_N),'MarkerSize',80,'LineWidth',1);  
    h(i).DisplayName = [datestr(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),'yyyy-mm-dd HH:MM') ...
                    ' | ' num2str(PcStructOut(i).PC2DBase,'%0.2E')];
end
hold off
box on
grid on
xlabel({'Crosstrack Change'; '(m)'})
ylabel({'Intrack Change';'(m)'})
Leg = legend(h,'LineWidth',1,'NumColumns',1);
title(Leg,'       Message Time  |     Pc   ');
legend(h);

saveas(gca,fullfile(FiguresFolder,'SecondaryRadialPlaneHistory.fig'));
saveas(gca,fullfile(FiguresFolder,'SecondaryRadialPlaneHistory.png'));
%% Generate Combined Conjunction RIC Uncertainty History
figure('NumberTitle','off','Name','RIC Miss/Uncertainty and Mahalanobis Update Plots','Position',[500,150,650,800])
% Mahalonobis Distance
subplot(4,1,1)
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Conjunction Mahalanobis Distance and RIC Miss Components Uncertainties'});
hold on
i=1;
% SecondaryDiff = zeros(length(CDMfile),3);
% SecondaryDiff(i,:) = zeros(3,1);
% PrimaryDiff = zeros(length(CDMfile),3);
% RICCombDiff = zeros(length(CDMfile),3);
% PrimaryDiff(i,:) = zeros(3,1);
MDCorr = zeros(length(CDMfile)-1,3);
clear h
for i = 2:length(CDMfile)
%     SecondaryDiff(i,:) = ([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z] + ...
%                           (datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60* ...
%                           [PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT]) - ...
%                           [PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z];
%     PrimaryDiff(i,:)   = ([PcStructOut(i).CDM.OBJECT1_X PcStructOut(i).CDM.OBJECT1_Y PcStructOut(i).CDM.OBJECT1_Z] + ...
%                           (datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60* ...
%                           [PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT]) - ...
%                           [PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z];
%     RICUpdateDiff(i,1:3) = ECI2RIC(SecondaryDiff(i,:),[PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z],[PcStructOut(i-1).CDM.OBJECT2_X_DOT PcStructOut(i-1).CDM.OBJECT2_Y_DOT PcStructOut(i-1).CDM.OBJECT2_Z_DOT])*1000;
%     RICUpdateDiff(i,4:6) = ECI2RIC(PrimaryDiff(i,:),[PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z],[PcStructOut(i-1).CDM.OBJECT1_X_DOT PcStructOut(i-1).CDM.OBJECT1_Y_DOT PcStructOut(i-1).CDM.OBJECT1_Z_DOT])*1000;
    MDCorr(i-1,1)  = sqrt(SecondaryDiff(i,:)*inv(PcStructOut(i-1).CDM.OBJECT2_C_J2K(1:3,1:3)/1000000)*SecondaryDiff(i,:)');
    MDCorr(i-1,2)  = sqrt(PrimaryDiff(i,:)*inv(PcStructOut(i-1).CDM.OBJECT1_C_J2K(1:3,1:3)/1000000)*PrimaryDiff(i,:)');
    RICCombDiff(i,:)    = [PcStructOut(i).CDM.RELATIVE_POSITION_R PcStructOut(i).CDM.RELATIVE_POSITION_T PcStructOut(i).CDM.RELATIVE_POSITION_N]...
                           - [PcStructOut(i-1).CDM.RELATIVE_POSITION_R PcStructOut(i-1).CDM.RELATIVE_POSITION_T PcStructOut(i-1).CDM.RELATIVE_POSITION_N];
    MDCorr(i-1,3)  = sqrt(RICCombDiff(i,:)*inv(PcStructOut(i-1).CDM.COMBINED_C_RTN(1:3,1:3))*RICCombDiff(i,:)');
    RICCombDiff(i,:)    = [PcStructOut(i).CDM.RELATIVE_POSITION_R PcStructOut(i).CDM.RELATIVE_POSITION_T PcStructOut(i).CDM.RELATIVE_POSITION_N];
    MDCorr(i-1,3)  = sqrt(RICCombDiff(i,:)*inv(PcStructOut(i).CDM.COMBINED_C_RTN(1:3,1:3))*RICCombDiff(i,:)');
    MDCorr(i-1,3)  = sqrt(RICCombDiff(i,:)*inv(diag(diag(PcStructOut(i).CDM.COMBINED_C_RTN(1:3,1:3))))*RICCombDiff(i,:)');
    if i==2
%         h(1) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,1),80,'rs','filled');  
%         h(2) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,2),80,'bo','filled');  
        h(1) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,3),80,'g^','filled');  
%         h(1).DisplayName = 'Secondary';
%         h(2).DisplayName = 'Primary';
        h(1).DisplayName = 'Combined';
    else
%         scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,1),80,'rs','filled');  
%         scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,2),80,'bo','filled');  
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,3),80,'g^','filled');  
    end
end
hold off
box on
grid on
ylabel({'Mahalanobis'; 'Distance'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
if max(MDCorr(:,3))>10
    set(gca, 'YScale', 'log')
    ylim([10^floor(log10(min(min(MDCorr(:,3))))) 10^ceil(log10(max(max(MDCorr(:,3)))))])
end

% legend(h);

% Radial Miss and Combined Uncertainty
subplot(4,1,2)
hold on
for i = 1:length(CDMfile)
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.RELATIVE_POSITION_R,sqrt(PcStructOut(i).CDM.COMBINED_C_RTN(1,1)),'s','Color',[0 0.4470 0.7410],'MarkerSize',5,'MarkerEdgeColor',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',1);
end
plot([floor(min(MessageTimes)) ceil(max(MessageTimes))],[0 0],'-k','LineWidth',0.5)
hold off
ylabel({'Radial Update'; 'Sequence (m)'})
box on
grid on
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')

% Intrack Miss and Combined Uncertainty
subplot(4,1,3)
hold on
for i = 1:length(CDMfile)
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.RELATIVE_POSITION_T,sqrt(PcStructOut(i).CDM.COMBINED_C_RTN(2,2)),'s','Color',[0.8500 0.3250 0.0980],'MarkerSize',5,'MarkerEdgeColor',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',1);
end
plot([floor(min(MessageTimes)) ceil(max(MessageTimes))],[0 0],'-k','LineWidth',0.5)
hold off
box on
grid on
ylabel({'Intrack Update'; 'Sequence (m)'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')

% Crosstrack miss and combined uncertainty
subplot(4,1,4)
hold on
for i = 1:length(CDMfile)
    errorbar(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.RELATIVE_POSITION_N,sqrt(PcStructOut(i).CDM.COMBINED_C_RTN(3,3)),'s','Color',[0.4660 0.6740 0.1880],'MarkerSize',5,'MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerFaceColor',[0.4660 0.6740 0.1880],'LineWidth',1);
end
plot([floor(min(MessageTimes)) ceil(max(MessageTimes))],[0 0],'-k','LineWidth',0.5)
hold off
box on
grid on
ylabel({'CrossTrack Update'; 'Sequence (m)'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
xlabel('Date of Message Update');

saveas(gca,fullfile(FiguresFolder,'CombinedRICUncertainty.fig'));
saveas(gca,fullfile(FiguresFolder,'CombinedRICUncertainty.png'));

%% Create PC/2D Mahalanobis Distance history Plot
figure('NumberTitle','off','Name','Pc and 2D Mahalanobis Plots','Position',[500,150,650,800])
% Mahalonobis Distance
subplot(2,1,1)
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Conjunction Plane Mahalanobis Distance and Pc Estimates'});
clear MDCorr
MDCorr = zeros(length(CDMfile),2);
for i = 1:length(CDMfile)
    MDCorr(i,1) = datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS');
    MDCorr(i,2) = sqrt([PcStructOut(i).CDM.MISS_DISTANCE 0]*inv(PcStructOut(i).CombinedCov_2D)*[PcStructOut(i).CDM.MISS_DISTANCE 0]');
end
scatter(MDCorr(:,1),MDCorr(:,2),80,'g^','filled');
box on
grid on
ylabel({'2D Mahalanobis'; 'Distance'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
if max(MDCorr(:,2))>10
    set(gca, 'YScale', 'log')
end
ylim([10^floor(log10(min(MDCorr(:,2)))) 10^ceil(log10(max(MDCorr(:,2))))])
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Conjunction Plane Mahalanobis Distance and Pc Estimates'});

% legend(h);

% Pc Estimates
subplot(2,1,2)
hold on
scatter(MDCorr(:,1),max([PcStructOut(:).PC3DBase],1E-10),80,'filled');
scatter(MDCorr(:,1),max([PcStructOut(:).DilutionMaxPcInputHBR],1E-10),80,'filled');
hold off
box on
grid on
ylabel({'Probability of Collision'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
xlabel('Date of Message Update');
legend('3D Pc','Dilution Max Pc','SouthWest');
set(gca, 'YScale', 'log')

saveas(gca,fullfile(FiguresFolder,'MD2D_Pc_History.fig'));
saveas(gca,fullfile(FiguresFolder,'MD2D_Pc_History.png'));
%% Generate Secondary Object EDR/WRMS Uncertainty History
figure('NumberTitle','off','Name','EDR, WRMS and Mahalanobis Update Plots','Position',[500,150,650,800])
% Mahalonobis Distance
subplot(3,1,1)
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Mahalanobis Distance, EDR, and WRMS Update Cycles'});
hold on
i=1;
% SecondaryDiff = zeros(length(CDMfile),3);
% SecondaryDiff(i,:) = zeros(3,1);
% PrimaryDiff = zeros(length(CDMfile),3);
% PrimaryDiff(i,:) = zeros(3,1);
MDCorr = zeros(length(CDMfile)-1,2);
for i = 2:length(CDMfile)
%     SecondaryDiff(i,:) = ([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z] + ...
%                           (datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) * 24*60*60*...
%                           [PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT]) - ...
%                           [PcStructOut(i-1).CDM.OBJECT2_X PcStructOut(i-1).CDM.OBJECT2_Y PcStructOut(i-1).CDM.OBJECT2_Z];
%     PrimaryDiff(i,:)   = ([PcStructOut(i).CDM.OBJECT1_X PcStructOut(i).CDM.OBJECT1_Y PcStructOut(i).CDM.OBJECT1_Z] + ...
%                           (datenum(PcStructOut(i-1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) * 24*60*60*...
%                           [PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT]) - ...
%                           [PcStructOut(i-1).CDM.OBJECT1_X PcStructOut(i-1).CDM.OBJECT1_Y PcStructOut(i-1).CDM.OBJECT1_Z];
    MDCorr(i-1,1)  = sqrt(SecondaryDiff(i,:)*inv(PcStructOut(i-1).CDM.OBJECT2_C_J2K(1:3,1:3)/1000000)*SecondaryDiff(i,:)');
    MDCorr(i-1,2)  = sqrt(PrimaryDiff(i,:)*inv(PcStructOut(i-1).CDM.OBJECT1_C_J2K(1:3,1:3)/1000000)*PrimaryDiff(i,:)');
    scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,1),80,'rs','filled');  
    scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),MDCorr(i-1,2),80,'bo','filled');  
end
hold off
box on
grid on
ylabel({'Mahalanobis'; 'Correction'})
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
set(gca, 'YScale', 'log')
ylim([10^floor(log10(min(min(MDCorr(:,1:2))))) 10^ceil(log10(max(max(MDCorr(:,1:2)))))])

subplot(3,1,2)
clear h
hold on
for i = 1:length(CDMfile)
    if i==1
        h(1) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_SEDR,80,'rs','filled');
        h(2) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_SEDR,80,'bo','filled');
        h(1).DisplayName = 'Secondary';
        h(2).DisplayName = 'Primary';
    else
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_SEDR,80,'rs','filled');
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_SEDR,80,'bo','filled');
    end
end
hold off
ylabel({'EDR'; 'W/kg'})
box on
grid on
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')
legend(h)

subplot(3,1,3)
hold on
for i = 1:length(CDMfile)
    scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_WEIGHTED_RMS,80,'rs','filled');
    scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_WEIGHTED_RMS,80,'bo','filled');
end
hold off
ylabel('WRMS')
box on
grid on
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd HH:MM')

saveas(gca,fullfile(FiguresFolder,'EDRHistory.fig'));
saveas(gca,fullfile(FiguresFolder,'EDRHistory.png'));
%% Generate Secondary Object Space Weather Sensitivity
figure('NumberTitle','off','Name','Space Weather Sensitivity Plots','Position',[500,150,650,800])
% Build Space Weather Entries
SpaceIdx = find([SWData(:).Datenum] >= floor(datenum(strrep(PcStructOut(1).CDM.CREATION_DATE,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')) & [SWData(:).Datenum] <= ceil(datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')));
SWPlotData = zeros((length(SpaceIdx)-1)*8+1,5);
SWPlotData(:,1) = SWData(SpaceIdx(1)).Datenum:1/8:SWData(SpaceIdx(end)).Datenum;
startIdx = 1;
for i=1:length(SpaceIdx)-1
    SWPlotData(startIdx:startIdx+7,2) = [SWData(SpaceIdx(i)).Ap0003
                                         SWData(SpaceIdx(i)).Ap0306
                                         SWData(SpaceIdx(i)).Ap0609
                                         SWData(SpaceIdx(i)).Ap0912
                                         SWData(SpaceIdx(i)).Ap1215
                                         SWData(SpaceIdx(i)).Ap1518
                                         SWData(SpaceIdx(i)).Ap1821
                                         SWData(SpaceIdx(i)).Ap2100];
    SWPlotData(startIdx:startIdx+7,3) = SWData(SpaceIdx(i)).ApAvg;
    SWPlotData(startIdx:startIdx+7,4) = SWData(SpaceIdx(i)).F107Obs;
    SWPlotData(startIdx:startIdx+7,5) = SWData(SpaceIdx(i)).F107Obs81Ctr;
    startIdx = startIdx+8;
end
SWPlotData(end,2:5) = SWPlotData(end-1,2:5);
clear startIdx

% Space Weather
subplot(4,1,1)
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Space Weather Sensitivity Plots'});
yyaxis left
hold on
plot(SWPlotData(:,1),SWPlotData(:,2),'LineWidth',2)
plot(SWPlotData(:,1),SWPlotData(:,3),'LineWidth',1)
hold off
ylabel({'Geomagnetic Activity'; 'Ap'})
ylim([0 max(SWPlotData(:,2))+5]);
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
yyaxis right
hold on
plot(SWPlotData(:,1),SWPlotData(:,4),'LineWidth',2)
plot(SWPlotData(:,1),SWPlotData(:,5),'LineWidth',1)
hold off
ylim([0 max(SWPlotData(:,4))+10]);
ylabel({'Solar Activity'; 'F_1_0_._7'})
box on
grid on
legend('3-Hourly Ap','Daily Avg Ap','Obs Flux','81Day Avg Flux')

% xlim([SWPlotData(1,1) SWPlotData(end,1)]);
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd')

% Secondary Radial Miss Update Cycle changes
subplot(4,1,2)
RICUpdateDiff = [zeros(size(SecondaryDiff)) zeros(size(SecondaryDiff))];
hold on
for i = 1:length(CDMfile)
    [r_s2,v_s2] = PosVelConvert([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z],...
                                [PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT],...
                                strrep(PcStructOut(i).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    [r_s1,v_s1] = PosVelConvert([PcStructOut(1).CDM.OBJECT2_X PcStructOut(1).CDM.OBJECT2_Y PcStructOut(1).CDM.OBJECT2_Z],...
                                [PcStructOut(1).CDM.OBJECT2_X_DOT PcStructOut(1).CDM.OBJECT2_Y_DOT PcStructOut(1).CDM.OBJECT2_Z_DOT],...
                                strrep(PcStructOut(1).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    SecondaryDiff(i,:) = (r_s2 ...
                          +(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
                          *v_s2 ...
                          )-r_s1;
    RICUpdateDiff(i,1:3) = ECI2RIC(SecondaryDiff(i,:),r_s1,v_s1)*1000;
    [r_s2,v_s2] = PosVelConvert([PcStructOut(i).CDM.OBJECT1_X PcStructOut(i).CDM.OBJECT1_Y PcStructOut(i).CDM.OBJECT1_Z],...
                                [PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT],...
                                strrep(PcStructOut(i).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    [r_s1,v_s1] = PosVelConvert([PcStructOut(1).CDM.OBJECT1_X PcStructOut(1).CDM.OBJECT1_Y PcStructOut(1).CDM.OBJECT1_Z],...
                                [PcStructOut(1).CDM.OBJECT1_X_DOT PcStructOut(1).CDM.OBJECT1_Y_DOT PcStructOut(1).CDM.OBJECT1_Z_DOT],...
                                strrep(PcStructOut(1).CDM.TCA,'T',' '),'ECF2J2K','4terms');
    PrimaryDiff(i,:)   = (r_s2 ...
                          +(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) *24*60*60 ...
                          *v_s2 ...
                          )-r_s1;
    RICUpdateDiff(i,4:6) = ECI2RIC(PrimaryDiff(i,:),r_s1,v_s1)*1000;
%     SecondaryDiff(i,:) = ([PcStructOut(i).CDM.OBJECT2_X PcStructOut(i).CDM.OBJECT2_Y PcStructOut(i).CDM.OBJECT2_Z] + ...
%                           (datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) * 24*60*60*...
%                           [PcStructOut(i).CDM.OBJECT2_X_DOT PcStructOut(i).CDM.OBJECT2_Y_DOT PcStructOut(i).CDM.OBJECT2_Z_DOT]) - ...
%                           [PcStructOut(1).CDM.OBJECT2_X PcStructOut(1).CDM.OBJECT2_Y PcStructOut(1).CDM.OBJECT2_Z];
%     PrimaryDiff(i,:) =  ([PcStructOut(i).CDM.OBJECT1_X PcStructOut(i).CDM.OBJECT1_Y PcStructOut(i).CDM.OBJECT1_Z] + ...
%                           (datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')-datenum(PcStructOut(i).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF')) * 24*60*60*...
%                           [PcStructOut(i).CDM.OBJECT1_X_DOT PcStructOut(i).CDM.OBJECT1_Y_DOT PcStructOut(i).CDM.OBJECT1_Z_DOT]) - ...
%                           [PcStructOut(1).CDM.OBJECT1_X PcStructOut(1).CDM.OBJECT1_Y PcStructOut(1).CDM.OBJECT1_Z];
%     RICUpdateDiff(i,1:3) = ECI2RIC(SecondaryDiff(i,:),[PcStructOut(1).CDM.OBJECT2_X PcStructOut(1).CDM.OBJECT2_Y PcStructOut(1).CDM.OBJECT2_Z],[PcStructOut(1).CDM.OBJECT2_X_DOT PcStructOut(1).CDM.OBJECT2_Y_DOT PcStructOut(1).CDM.OBJECT2_Z_DOT])*1000;
%     RICUpdateDiff(i,4:6) = ECI2RIC(PrimaryDiff(i,:),[PcStructOut(1).CDM.OBJECT1_X PcStructOut(1).CDM.OBJECT1_Y PcStructOut(1).CDM.OBJECT1_Z],[PcStructOut(1).CDM.OBJECT1_X_DOT PcStructOut(1).CDM.OBJECT1_Y_DOT PcStructOut(1).CDM.OBJECT1_Z_DOT])*1000;
    scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,1),80,'rs','filled');  
    scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),RICUpdateDiff(i,4),80,'bo','filled');  
end
plot([SWPlotData(1,1) SWPlotData(end,1)],[0 0],'-k','LineWidth',0.5)
hold off
ylabel({'Secondary';'Radial Difference'; 'from Initial (m)'})
box on
grid on
% xlim([SWPlotData(1,1) SWPlotData(end,1)]);
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))]);
datetick('x','mm-dd')
% datetick('x','mm-dd HH:MM')

% Plot EDR and BC for the secondary
subplot(4,1,3)
yyaxis left
clear h
hold on
for i = 1:length(CDMfile)
    if i==1
        h(1) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_SEDR,80,'s','filled');
        h(2) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_SEDR,80,'o','filled');
        h(1).DisplayName = 'Secondary EDR';
        h(2).DisplayName = 'Primary EDR';
    else
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_SEDR,80,'s','filled');
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_SEDR,80,'o','filled');
    end
end
hold off
ylabel({'EDR'; 'W/kg'})
% set(gca, 'YScale', 'log')
yyaxis right
hold on
for i = 1:length(CDMfile)
    if i==1
        h(3) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_CD_AREA_OVER_MASS,80,'s','filled');
        h(4) = scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_CD_AREA_OVER_MASS,80,'o','filled');
        h(3).DisplayName = 'Secondary BC';
        h(4).DisplayName = 'Primary BC';
    else
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_CD_AREA_OVER_MASS,80,'s','filled');
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_CD_AREA_OVER_MASS,80,'o','filled');
    end
end
hold off
ylabel({'BC'; 'm^2/kg'})
set(gca, 'YScale', 'log')
box on
grid on
% xlim([SWPlotData(1,1) SWPlotData(end,1)]);
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))+3]);
datetick('x','mm-dd')
% datetick('x','mm-dd HH:MM')
Leg = legend(h,'LineWidth',1);

% Plot DCP Parameters
subplot(4,1,4)
hold on
for i = 1:length(CDMfile)
    if isfield(PcStructOut(i).CDM,'OBJECT2_DCP_DENSITY_UNCERTAINTY') && ~isempty(PcStructOut(i).CDM.OBJECT2_DCP_DENSITY_UNCERTAINTY)
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT2_DCP_DENSITY_UNCERTAINTY*100,80,'rs','filled');
    else
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),0,80,'rs','filled');
    end
    if isfield(PcStructOut(i).CDM,'OBJECT1_DCP_DENSITY_UNCERTAINTY') && ~isempty(PcStructOut(i).CDM.OBJECT1_DCP_DENSITY_UNCERTAINTY)
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),PcStructOut(i).CDM.OBJECT1_DCP_DENSITY_UNCERTAINTY*100,80,'bo','filled');
    else
        scatter(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),0*100,80,'bo','filled');
    end
end
hold off
if isfield(PcStructOut(i).CDM,'OBJECT2_DCP_DENSITY_UNCERTAINTY')
    ylabel({'DCP \sigma';'%'})
else
    ylabel({'DCP \sigma';'(Not Available)'})
    set(gca,'Color',[0.75 0.75 0.75])
end
box on
grid on
legend('Secondary','Primary')
% xlim([SWPlotData(1,1) SWPlotData(end,1)]);
xlim([floor(min(MessageTimes)) ceil(max(MessageTimes))+3]);
datetick('x','mm-dd')
% datetick('x','mm-dd HH:MM')

saveas(gca,fullfile(FiguresFolder,'SWSensitivity.fig'));
saveas(gca,fullfile(FiguresFolder,'SWSensitivity.png'));

%% Generate Conjunction Plane Plots
% Creating figure
figure('NumberTitle','off','Name','2D Conjunction Plane','Position',[500,250,840,700]);
clear h
hold on
for i = 1:length(CDMfile)
    ConjData = [PcStructOut(i).SemiMajor2D ...
                PcStructOut(i).SemiMinor2D ...
                PcStructOut(i).ClockAngle ...
                PcStructOut(i).InputHBR/1000 ...
                PcStructOut(i).CDM.MISS_DISTANCE ...
                sqrt(PcStructOut(i).CombinedCov_2D(1,1))];
%     GenConjPlanePlot(PcStructOut(i).PC2DBase,ConjData)
    a          = PcStructOut(i).SemiMajor2D;
    b          = PcStructOut(i).SemiMinor2D;
    ClockAngle = PcStructOut(i).ClockAngle;
    HBR        = PcStructOut(i).InputHBR/1000;
    MissDist   = PcStructOut(i).CDM.MISS_DISTANCE;
    MinHBSigmaDistance = (MissDist - HBR) / ConjData(6);
    
    % Generate basis ellipse (parameterization)
    u  = (0:0.1:2*pi+0.1)';
    x0 = a * cos(u);
    y0 = b * sin(u);
    
    % Rotate ellipse by clock angle
    x = x0.*cosd(ClockAngle) - y0.*sind(ClockAngle); 
    y = x0.*sind(ClockAngle) + y0.*cosd(ClockAngle);
    
    % Construct HBR
    xHBR = HBR * cos(u) + MissDist;
    yHBR = HBR * sin(u);
    
    % Plot 1, sqrt(2), and 3-sigma ellipses
    h(i) = plot(sqrt(2)*x, sqrt(2)*y,'LineWidth',1);
%     h(i).DisplayName = ['Message Time: ' datestr(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),'yyyy-mm-dd HH:MM') ...
%                      ' | Pc = ' num2str(PcStructOut(i).PC2DBase,'%0.2E')];
    h(i).DisplayName = [datestr(datenum(PcStructOut(i).CDM.CREATION_DATE,'yyyy-mm-ddTHH:MM:SS'),'yyyy-mm-dd HH:MM') ...
                        ' | ' num2str(PcStructOut(i).PC2DBase,'%0.2E')];
    
    % Plot HBR
%     plot(xHBR,yHBR,'Color',get(h,'Color'),'LineWidth',2);
    scatter(MissDist,0,[],get(h(i),'Color'),'filled')
    
    
end
hold off
LegendHBR = strcat('Hard Body Region (',num2str(HBR*1000),' m)');
title({['Conjunction: ' PrimaryText ' vs ' PcStructOut(1).CDM.OBJECT2_OBJECT_DESIGNATOR ' | TCA: ' datestr(datenum(PcStructOut(1).CDM.TCA,'yyyy-mm-ddTHH:MM:SS.FFF'),'yyyy-mm-dd HH:MM')]; 'Conjunction Plane - SQRT Two Sigma Covariance'});
xlabel('Miss Distance Direction [km]');
ylabel('Relative Out of Plane [km]');
% legend('One Sigma Covariance', 'SQRT Two Sigma Covariance', 'Three Sigma Covariance',LegendHBR,'Location','NorthEast');
Leg = legend(h,'LineWidth',1);
title(Leg,'       Message Time  |     Pc   ')
box on
grid on

% Correct Scale of Y axis
if max(abs(get(gca,'ylim')))<min(abs(get(gca,'xlim')))
    % Set y scale to be comparable to x scale (in negative direction)
    ylim([-min(abs(get(gca,'xlim'))) min(abs(get(gca,'xlim')))]);
else
    % increase scale of xlimits
    xlim(get(gca,'ylim'));
end

% save figures
saveas(gca,fullfile(FiguresFolder,'ConjunctionPlane.fig'));
saveas(gca,fullfile(FiguresFolder,'ConjunctionPlane.png'));

%% Generate Space Weather Prediction Plots
clear h
if ~isempty(SolarFluxfile)
    % Set Lower limit on dates to plot
    DateLimit = floor(MessageTimes(1))-3;

    % Parse all input Solar_flux files
    for ii = 1:length(SolarFluxfile)
        [JBHSGI_Data{ii},~] = Parse_jbhsgi(SolarFluxfile{ii});
    end

    figure('NumberTitle','off','Name','JBHSGI Solar Data','Position',[500,150,700,840]);
    subplot(4,1,1)
    hold on
        % Plot F10.7
        for ii = 1:length(JBHSGI_Data)
            idx1 = find(JBHSGI_Data{ii}(:,13)==0);
            idx2 = find(JBHSGI_Data{ii}(:,13)==1); idx2 = [idx2(1)-1; idx2];
            idx3 = find(JBHSGI_Data{ii}(:,13)==2); idx3 = [idx3(1)-1; idx3];
            h(ii) = plot(JBHSGI_Data{ii}(idx1,1),JBHSGI_Data{ii}(idx1,2),'-' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx2,1),JBHSGI_Data{ii}(idx2,2),'-.','Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx3,1),JBHSGI_Data{ii}(idx3,2),':' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
            CreationTimeIdx = find(JBHSGI_Data{ii}(:,13)==2,1,'first');
            h(ii).DisplayName = datestr(JBHSGI_Data{ii}(CreationTimeIdx,1),'yyyy-mm-dd HH:MM');
        end
    hold off
    xlim([DateLimit ceil(datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')+0.5)]);
    yl = ylim;
    ylim([0 yl(2)])
    ylabel('F_1_0_._7')
    datetick('x','mmmdd','keeplimits')
    box on
    grid on
    title({'JBHSGI Solar Flux Data File Predictions and Estimates'; ...
            'Solid Lines - Historical Values'; ...
            'Dash-Dotted Lines - Estimated Values'; ...
            'Dotted Lines - Predicted Values'})

    subplot(4,1,2)
    hold on
        % Plot Ap
        for ii = 1:length(JBHSGI_Data)
            idx1 = find(JBHSGI_Data{ii}(:,13)==0);
            idx2 = find(JBHSGI_Data{ii}(:,13)==1); idx2 = [idx2(1)-1; idx2];
            idx3 = find(JBHSGI_Data{ii}(:,13)==2); idx3 = [idx3(1)-1; idx3];
            h(ii) = plot(JBHSGI_Data{ii}(idx1,1),JBHSGI_Data{ii}(idx1,10),'-' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx2,1),JBHSGI_Data{ii}(idx2,10),'-.','Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx3,1),JBHSGI_Data{ii}(idx3,10),':' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
            CreationTimeIdx = find(JBHSGI_Data{ii}(:,13)==2,1,'first');
            h(ii).DisplayName = datestr(JBHSGI_Data{ii}(CreationTimeIdx,1),'yyyy-mm-dd HH:MM');
        end
    hold off
    xlim([DateLimit ceil(datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')+0.5)]);
    yl = ylim;
    ylim([0 yl(2)])
    Leg = legend(h,'LineWidth',1);
    ylabel('Ap')
    datetick('x','mmmdd','keeplimits')
    box on
    grid on

    subplot(4,1,3)
    hold on
        % Plot Dst
        for ii = 1:length(JBHSGI_Data)
            idx1 = find(JBHSGI_Data{ii}(:,13)==0);
            idx2 = find(JBHSGI_Data{ii}(:,13)==1); idx2 = [idx2(1)-1; idx2];
            idx3 = find(JBHSGI_Data{ii}(:,13)==2); idx3 = [idx3(1)-1; idx3];
            h(ii) = plot(JBHSGI_Data{ii}(idx1,1),JBHSGI_Data{ii}(idx1,11),'-' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx2,1),JBHSGI_Data{ii}(idx2,11),'-.','Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx3,1),JBHSGI_Data{ii}(idx3,11),':' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
            CreationTimeIdx = find(JBHSGI_Data{ii}(:,13)==2,1,'first');
            h(ii).DisplayName = datestr(JBHSGI_Data{ii}(CreationTimeIdx,1),'yyyy-mm-dd HH:MM');
        end
    hold off
    xlim([DateLimit ceil(datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')+0.5)]);
    yl = ylim;
    ylim([yl(1) 0])
    ylabel('Dst')
    datetick('x','mmmdd','keeplimits')
    box on
    grid on

    subplot(4,1,4)
    hold on
        % Plot Dtc
        for ii = 1:length(JBHSGI_Data)
            idx1 = find(JBHSGI_Data{ii}(:,13)==0);
            idx2 = find(JBHSGI_Data{ii}(:,13)==1); idx2 = [idx2(1)-1; idx2];
            idx3 = find(JBHSGI_Data{ii}(:,13)==2); idx3 = [idx3(1)-1; idx3];
            h(ii) = plot(JBHSGI_Data{ii}(idx1,1),JBHSGI_Data{ii}(idx1,12),'-' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx2,1),JBHSGI_Data{ii}(idx2,12),'-.','Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
                    plot(JBHSGI_Data{ii}(idx3,1),JBHSGI_Data{ii}(idx3,12),':' ,'Color',RGBTriplets(mod(ii-1,size(RGBTriplets,1))+1,:),'LineWidth',1.5);
            CreationTimeIdx = find(JBHSGI_Data{ii}(:,13)==2,1,'first');
            h(ii).DisplayName = datestr(JBHSGI_Data{ii}(CreationTimeIdx,1),'yyyy-mm-dd HH:MM');
        end
    hold off
    xlim([DateLimit ceil(datenum(strrep(PcStructOut(1).CDM.TCA,' ','T'),'yyyy-mm-ddTHH:MM:SS.FFF')+0.5)]);
    yl = ylim;
    ylim([0 yl(2)])
    ylabel('dtc')
    datetick('x','mmmdd','keeplimits')
    box on
    grid on


    saveas(gca,fullfile(FiguresFolder,'SpaceWeatherPrediction.fig'));
    saveas(gca,fullfile(FiguresFolder,'SpaceWeatherPrediction.png'));
end

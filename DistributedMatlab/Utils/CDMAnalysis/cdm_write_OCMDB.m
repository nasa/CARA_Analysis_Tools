function [success] = cdm_write_OCMDB(FolderName,CDMFilename,DB,DB_index,params)
% =========================================================================
%
% Function to write conjunction data message from an OCMDB row format 
%
% =========================================================================
%
% INPUT:
%
%   FolderName          = Folder to which to write the resultant CDM file
%   CDMFilename         = Output filename to which to write the resultant
%                         CDM File (if '' or [] is chosen as input, the
%                         program will automatically generate the conjID
%                         and use it as the CDM filename)
%   DB                  = An OCMDB data matrix           [nx218] or [nx234]
%   DB_index            = index (row number) of event to be written to CDM 
%                         text file (the conj object may contain several 
%                         different events, hence, a specific event must be 
%                         identified for writing) (optional - default is
%                         set to 1)
%   params              = Structure containing optional inputs. The
%                         possible params fields are provided in the
%                         following. (optional)
%
% OUTPUT:
%
%   success             = Boolean output indicating success(1) or failure
%                         (0) of the CDM writing process
%
% -------------------------------------------------------------------------
%
%   Fields within the optional params structure:
%
%   params.priName      = Name of the primary object if available
%                         (Default primary name is set to "UNKNOWN")
%                         the generated CDM)
%   params.secName      = Name of the secondary object if available
%                         (Default secondary name is set to "UNKNOWN")
%   params.priIntDes    = International designator of the primary
%                         (Default set to "UNKNOWN")
%   params.secIntDes    = International designator of the secondary
%                         (Default set to "UNKNOWN")
%   params.ref_mode     = Reference frame for the object states provided
%                         in the CDM. 1 - EME2000, 2 - ITRF
%                         (Default reference frame is set to 2)
%   params.ooInfo.priIsOO
%                       = Flag representing whether primary data are
%                         provided through Operator/Owner ephemerides 
%                         (Default set to false)
%   params.ooInfo.secIsOO
%                       = Flag representing whether secondary data are
%                         provided through Operator/Owner ephemerides 
%                         (Default set to false)
%
% =========================================================================


% Set paths
persistent pathsLoaded
if isempty(pathsLoaded)
    [mpath,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(mpath,'../PosVelTransformations')); addpath(s.path);
    s = what(fullfile(mpath,'../General')); addpath(s.path);
    pathsLoaded = true;
end

% Check number of inputs
Nargin = nargin;
if Nargin == 3
    DB_index = 1;
    params = [];
elseif Nargin == 4
    params = [];
elseif Nargin ~= 5
    error('Incorrect number of inputs')
end

% Setting default parameters
params = set_default_param(params,'priName','UNKNOWN');
params = set_default_param(params,'secName','UNKNOWN');
params = set_default_param(params,'priIntDes','UNKNOWN');
params = set_default_param(params,'secIntDes','UNKNOWN');
params = set_default_param(params,'ref_mode',2);
params = set_default_param(params,'ooInfo',[]);
params.ooInfo =  set_default_param(params.ooInfo,'priIsOO',false);
params.ooInfo =  set_default_param(params.ooInfo,'secIsOO',false);


dateFmt = 'yyyyMMdd_HHmmss';
% Loop through and create CDMs for each conjunction in the DB
i  = DB_index;
% Generate the conjunction ID
primaryID = DB(i,1);
secondaryID = DB(i,2);
primaryStr = num2str(primaryID,'%09d');
secondaryStr = num2str(secondaryID,'%09d');
createTime = datetime(DB(i,3),DB(i,4),DB(i,5),DB(i,6),DB(i,7),DB(i,8));
tcaTime = datetime(DB(i,12),DB(i,13),DB(i,14),DB(i,15),DB(i,16),DB(i,17));
createTimeStr = char(createTime,dateFmt);
tcaTimeStr = char(tcaTime,dateFmt);
createTimeStr2 = GetCDMDateStr(DB(i,3),DB(i,4),DB(i,5),DB(i,6),DB(i,7),DB(i,8),DB(i,9));
tcaTimeStr2 = GetCDMDateStr(DB(i,12),DB(i,13),DB(i,14),DB(i,15),DB(i,16),DB(i,17),DB(i,18));
conjID = [primaryStr '_conj_' secondaryStr '_' tcaTimeStr '_' createTimeStr];

% Get the object names and international designators for the
% primary and secondary
primaryName = params.priName;
secondaryName = params.secName;
primaryIntlDes = params.priIntDes;
secondaryIntlDes = params.secIntDes;

% Get the O/O indicator
if params.ooInfo.priIsOO || DB(i,37) <= 1969
    priODSource = 'O/O';
else
    priODSource = 'ASW';
end
if params.ooInfo.secIsOO || DB(i,97) <= 1969
    secODSource = 'O/O';
else
    secODSource = 'ASW';
end


% CDM file name
if isempty(CDMFilename)
    CDMFile = fullfile(FolderName, [conjID '.cdm']);
else
    file_name_parts = strsplit(CDMFilename,'.');
    CDMFile = fullfile(FolderName, [file_name_parts{1} '.cdm']);
end
fh = fopen(CDMFile,'w');
try
    % CDM Header Fields
    PrintCDMField(fh,'CCSDS_CDM_VERS','1.0','%s');
    PrintCDMField(fh,'CREATION_DATE',createTimeStr2,'%s');
    PrintCDMField(fh,'ORIGINATOR','CARA','%s'); % Previously JSPOC
    PrintCDMField(fh,'MESSAGE_FOR',primaryName,'%s');
    PrintCDMField(fh,'MESSAGE_ID',conjID,'%s');
    PrintCDMField(fh,'COMMENT SCREENING_OPTION','Covariance','%s');
    PrintCDMField(fh,'TCA',tcaTimeStr2,'%s');
    PrintCDMField(fh,'MISS_DISTANCE',DB(i,19),'%d','[m]');
    PrintCDMField(fh,'RELATIVE_SPEED',DB(i,20),'%d','[m/s]');
    PrintCDMField(fh,'RELATIVE_POSITION_R',DB(i,21),'%.1f','[m]');
    PrintCDMField(fh,'RELATIVE_POSITION_T',DB(i,22),'%.1f','[m]');
    PrintCDMField(fh,'RELATIVE_POSITION_N',DB(i,23),'%.1f','[m]');
    PrintCDMField(fh,'RELATIVE_VELOCITY_R',DB(i,24),'%.1f','[m/s]');
    PrintCDMField(fh,'RELATIVE_VELOCITY_T',DB(i,25),'%.1f','[m/s]');
    PrintCDMField(fh,'RELATIVE_VELOCITY_N',DB(i,26),'%.1f','[m/s]');
    PrintCDMField(fh,'COLLISION_PROBABILITY',DB(i,156),'%.3e');
    PrintCDMField(fh,'COLLISION_PROBABILITY_METHOD','FOSTER-1992','%s');
    PrintCDMField(fh,'COMMENT HBR',DB(i,157),'%.1f','[m]');

    % Object specific data
    for j = 1:2
        objectTxt = ['OBJECT' num2str(j)];
        if j == 1
            objectID = primaryStr;
            objectName = primaryName;
            objectIntlDes = primaryIntlDes;
            objectODSource = priODSource;
            offset = 0;
            dcpOffset = 0;
            ecioffset = 0;
        else
            objectID = secondaryStr;
            objectName = secondaryName;
            objectIntlDes = secondaryIntlDes;
            objectODSource = secODSource;
            offset = 60;
            dcpOffset = 8;
            ecioffset = 6;
        end
        PrintCDMField(fh,'OBJECT',objectTxt,'%s');
        PrintCDMField(fh,'OBJECT_DESIGNATOR',objectID,'%s');
        PrintCDMField(fh,'CATALOG_NAME','SATCAT','%s');
        PrintCDMField(fh,'OBJECT_NAME',objectName,'%s');
        PrintCDMField(fh,'INTERNATIONAL_DESIGNATOR',objectIntlDes,'%s');
        PrintCDMField(fh,'EPHEMERIS_NAME','NONE','%s');
        PrintCDMField(fh,'COVARIANCE_METHOD','CALCULATED','%s');
        PrintCDMField(fh,'MANEUVERABLE','N/A','%s');
        if params.ref_mode == 2
            PrintCDMField(fh,'REF_FRAME','ITRF','%s');
        else
            PrintCDMField(fh,'REF_FRAME','EME2000','%s');
        end
        PrintCDMField(fh,'GRAVITY_MODEL',ocmGeop2ZTnumSTR(DB(i,59+offset)),'%s');
        PrintCDMField(fh,'ATMOSPHERIC_MODEL',GetAtmoModel(DB(i,60+offset)),'%s');
        PrintCDMField(fh,'N_BODY_PERTURBATIONS',GetNBody(DB(i,61+offset)),'%s');
        PrintCDMField(fh,'SOLAR_RAD_PRESSURE',GetYesNo(DB(i,62+offset)),'%s');
        PrintCDMField(fh,'EARTH_TIDES',GetYesNo(DB(i,63+offset)),'%s');
        PrintCDMField(fh,'INTRACK_THRUST',GetYesNo(DB(i,64+offset)),'%s');
        PrintCDMField(fh,'COMMENT COVARIANCE_SCALE_FACTOR',DB(i,65+offset),'%.3f');
        PrintCDMField(fh,'COMMENT EXCLUSION_VOLUME_RADIUS',DB(i,66+offset),'%d','[m]');
        if ~isnan(DB(i,37+offset))
            lastObDateStr = GetCDMDateStr(DB(i,37+offset),DB(i,38+offset),DB(i,39+offset),...
                DB(i,40+offset),DB(i,41+offset),DB(i,42+offset),DB(i,43+offset));
            PrintCDMField(fh,'TIME_LASTOB_START',lastObDateStr,'%s');
            PrintCDMField(fh,'TIME_LASTOB_END',lastObDateStr,'%s');
        end
        PrintCDMField(fh,'COMMENT OD_DATA_SOURCE',objectODSource,'%s');
        PrintCDMField(fh,'RECOMMENDED_OD_SPAN',DB(i,44+offset),'%.2f','[d]');
        PrintCDMField(fh,'ACTUAL_OD_SPAN',DB(i,45+offset),'%.2f','[d]');
        PrintCDMField(fh,'OBS_AVAILABLE',DB(i,47+offset),'%d');
        PrintCDMField(fh,'OBS_USED',DB(i,48+offset),'%d');
        PrintCDMField(fh,'TRACKS_AVAILABLE',DB(i,49+offset),'%d');
        PrintCDMField(fh,'TRACKS_USED',DB(i,50+offset),'%d');
        PrintCDMField(fh,'RESIDUALS_ACCEPTED',DB(i,46+offset),'%.1f','[%%]');
        PrintCDMField(fh,'WEIGHTED_RMS',DB(i,55+offset),'%.3f');
        PrintCDMField(fh,'COMMENT Apogee Altitude',DB(i,51+offset),'%d','[km]');
        PrintCDMField(fh,'COMMENT Perigee Altitude',DB(i,52+offset),'%d','[km]');
        PrintCDMField(fh,'COMMENT Inclination',DB(i,53+offset),'%.1f','[deg]');
        if ~isnan(DB(i,54+offset))
            PrintCDMField(fh,'AREA_PC',DB(i,54+offset),'%.4f','[m**2]');
        end
        PrintCDMField(fh,'CD_AREA_OVER_MASS',DB(i,56+offset),'%.6f','[m**2/kg]');
        PrintCDMField(fh,'CR_AREA_OVER_MASS',DB(i,57+offset),'%.6f','[m**2/kg]');
        PrintCDMField(fh,'SEDR',DB(i,58+offset),'%.6f','[W/kg]');
        if size(DB,2) > 218
            PrintCDMField(fh,'COMMENT DCP Density Forecast Uncertainty',DB(i,219+dcpOffset),'%.18e');
            dcpPosVecStr = sprintf('%.18e %.18e %.18e [m]',DB(i,221+dcpOffset:223+dcpOffset).*1000);
            dcpVelVecStr = sprintf('%.18e %.18e %.18e [m/sec]',DB(i,224+dcpOffset:226+dcpOffset).*1000);
            PrintCDMField(fh,'COMMENT DCP Sensitivity Vector RTN Pos',dcpPosVecStr,'%s');
            PrintCDMField(fh,'COMMENT DCP Sensitivity Vector RTN Vel',dcpVelVecStr,'%s');
        end
        EpochUTC = strrep(tcaTimeStr2,'T',' ');
        if params.ref_mode == 2
            [rITRF,vITRF] = PosVelConvert(DB(i,67+offset:69+offset)./1000,DB(i,70+offset:72+offset)./1000,EpochUTC,'TDR2ECF','4terms');
            PrintCDMField(fh,'X',rITRF(1),'%.18e','[km]');
            PrintCDMField(fh,'Y',rITRF(2),'%.18e','[km]');
            PrintCDMField(fh,'Z',rITRF(3),'%.18e','[km]');
            PrintCDMField(fh,'X_DOT',vITRF(1),'%.18e','[km/s]');
            PrintCDMField(fh,'Y_DOT',vITRF(2),'%.18e','[km/s]');
            PrintCDMField(fh,'Z_DOT',vITRF(3),'%.18e','[km/s]');
        else
            PrintCDMField(fh,'X', DB(i,172 + ecioffset) ,'%.18e','[km]');
            PrintCDMField(fh,'Y', DB(i,173 + ecioffset) ,'%.18e','[km]');
            PrintCDMField(fh,'Z', DB(i,174 + ecioffset) ,'%.18e','[km]');
            PrintCDMField(fh,'X_DOT', DB(i,175 + ecioffset) ,'%.18e','[km/s]');
            PrintCDMField(fh,'Y_DOT', DB(i,176 + ecioffset) ,'%.18e','[km/s]');
            PrintCDMField(fh,'Z_DOT', DB(i,177 + ecioffset) ,'%.18e','[km/s]');
        end
        PrintCDMField(fh,'CR_R',      DB(i,73+offset),'%.18e','[m**2]');
        PrintCDMField(fh,'CT_R',      DB(i,74+offset),'%.18e','[m**2]');
        PrintCDMField(fh,'CT_T',      DB(i,79+offset),'%.18e','[m**2]');
        PrintCDMField(fh,'CN_R',      DB(i,75+offset),'%.18e','[m**2]');
        PrintCDMField(fh,'CN_T',      DB(i,80+offset),'%.18e','[m**2]');
        PrintCDMField(fh,'CN_N',      DB(i,84+offset),'%.18e','[m**2]');
        PrintCDMField(fh,'CRDOT_R',   DB(i,76+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CRDOT_T',   DB(i,81+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CRDOT_N',   DB(i,85+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CRDOT_RDOT',DB(i,88+offset),'%.18e','[m**2/s**2]');
        PrintCDMField(fh,'CTDOT_R',   DB(i,77+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CTDOT_T',   DB(i,82+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CTDOT_N',   DB(i,86+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CTDOT_RDOT',DB(i,89+offset),'%.18e','[m**2/s**2]');
        PrintCDMField(fh,'CTDOT_TDOT',DB(i,91+offset),'%.18e','[m**2/s**2]');
        PrintCDMField(fh,'CNDOT_R',   DB(i,78+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CNDOT_T',   DB(i,83+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CNDOT_N',   DB(i,87+offset),'%.18e','[m**2/s]');
        PrintCDMField(fh,'CNDOT_RDOT',DB(i,90+offset),'%.18e','[m**2/s**2]');
        PrintCDMField(fh,'CNDOT_TDOT',DB(i,92+offset),'%.18e','[m**2/s**2]');
        PrintCDMField(fh,'CNDOT_NDOT',DB(i,93+offset),'%.18e','[m**2/s**2]');
    end
    success = true;
catch ME
    success = false;
    fclose(fh);
    disp('Failed to create CDM');
    
end
fclose(fh);


end

%% Additional Functions
function PrintCDMField(fh, fieldName, value, format, units)
    if nargin == 4
        units = '';
    elseif nargin ~= 5
        error('Incorrect number of arguments passed in');
    end
    
    % Set the field name format
    defaultFieldLength = 43;
    if length(fieldName) > defaultFieldLength
        fieldNameFormat = ['%-' num2str(length(fieldName)) 's'];
    else
        fieldNameFormat = ['%-' num2str(defaultFieldLength) 's'];
    end
    
    % Adjust the field name format for COMMENT lines
    if startsWith(fieldName,'COMMENT ')
        fieldNameFormat = '%s';
    end
    
    % Set the units text
    if ~isempty(units)
        unitsTxt = [' ' units];
    else
        unitsTxt = '';
    end
    
    % Convert to integer value if integer outputs are needed
    if endsWith(format,'d')
        value = int64(value);
    end
    
    % Print the line
    fprintf(fh,[fieldNameFormat ' = ' format unitsTxt '\n'],fieldName,value);
end


function Fullstr = ocmGeop2ZTnumSTR (GeopNum)
% Convert the Gravity Model description within the OCMDB row to minimum degree
% term in Zonal and Terrestrial model 

if GeopNum == 1 || GeopNum == 9 || GeopNum == 17
    minDegree = 0;
elseif GeopNum == 10
    minDegree = 7;
elseif GeopNum == 2 || GeopNum == 11
    minDegree = 8;
elseif GeopNum == 3 || GeopNum == 12
    minDegree = 12;
elseif GeopNum == 4 || GeopNum == 13
    minDegree = 18;
elseif GeopNum == 5 || GeopNum == 14
    minDegree = 24;
elseif GeopNum == 6 || GeopNum == 15
    minDegree = 36;
elseif GeopNum == 7 || GeopNum == 16
    minDegree = 48;
elseif GeopNum == 8
    minDegree = 70;
else
    minDegree = nan;
end

minDegree = num2str(minDegree);

if GeopNum == 0
    Fullstr = 'UNKNOWN: 0D 0O';
elseif GeopNum < 9
    Fullstr = ['EGM-96: ' minDegree 'D ' minDegree 'O'];
else
    Fullstr = ['Custom: ' minDegree 'D ' minDegree 'O'];
end

end

function [atmoModel] = GetAtmoModel(modelNum)
    if modelNum == 0
        atmoModel = 'NONE';
    elseif modelNum == 1
        atmoModel = 'JACCHIA70DCA';
    elseif modelNum == 2
        atmoModel = 'JBH09';
    else
        atmoModel = 'UNKNOWN';
    end
end

function [nBody] = GetNBody(modelNum)
    if modelNum == 0
        nBody = 'NONE';
    elseif modelNum == 1
        nBody = 'MOON,SUN';
    else
        nBody = 'UNKNOWN';
    end
end

function [yesNo] = GetYesNo(modelNum)
    if modelNum == 0
        yesNo = 'NO';
    elseif modelNum == 1
        yesNo = 'YES';
    else
        yesNo = 'NO';
    end
end

function [dateStr] = GetCDMDateStr(year, month, day, hour, min, sec, millisec)
    sec = sec + millisec / 1000;
    dateStr = sprintf('%4d-%02d-%02dT%02d:%02d:%06.3f', year, month, day, hour, min, sec);
end
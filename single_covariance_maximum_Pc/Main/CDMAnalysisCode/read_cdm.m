function [cdmhead, cdmobj, Status] = read_cdm(filename)
% Reads CCSDS CDM file
%
% Purpose:  Reads a CCSDS Conjunction Data Message (CDM) returning two
% data structures.  The first contains the CDM header and relative metadata
% information.  The second contains the object 1 and 2 data.  
%
% Invocation Method: [cdmhead, cdmobj] = read_cdm(filename)
%
% Argument       I/O  Description
% --------       ---  ----------------------------------------------------
% filename        i   Fully qualified name of CDM file to read
% cdmhead         o   Structure containing the CDM header and relative meta
%                     data
% cdmobj(2)       o   Structure containing the object data for each of the
%                     two objects
%
% External References:
% Function Name      Purpose
% -------------      ------------------------------------------------------
% isempty_cell       Returns a logical array indicating which cells in a
%                    cell array are empty
% getCcsdsTimeFormat Determines format of CCSDS time string
% DOY2Date           Converts day of year to MATLAB date number
%
% Internal Functions:
% Function Name     Purpose
% -------------     ------------------------------------------------------
% get_keyword_value Gets next keyword-value pair from the CCSDS file
%
%
% Global References:
% Parameters        Description
% -------------     ------------------------------------------------------
% None
%
% Copyright C {2016} United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% This software developed under RIGHTS IN DATA - special works
% (FAR 52.227-17) as modified by NFS 1852.227-17.

% Development History:
% Name          Date        Description of Change
% ------------  ----------  ----------------------------------------------
% R. Coon       06/14/2016  Created
% D. Hall       09/16/2019  Modified to parse special COMMENT fields
%                           for the recently revised CDM format.
% D. Hall       10/07/2019  Modified to parse additional special COMMENTs
%                           for the recently revised CDM format.

% Initialize output
cdmhead = [];
cdmobj  = [];
Status = 0;

% List of CDM header and relative meta data keywords
% CCSDS CDM Version keyword must be the first entry
% Column 1 is the keyword.  
% Column 2 is the value format.  's' - string, 'f' -  numeric, 't' - time
headParms = {
   'CCSDS_CDM_VERS', 's';
   'CREATION_DATE', 't';
   'ORIGINATOR', 's';
   'MESSAGE_FOR', 's';
   'MESSAGE_ID', 's';
   'TCA', 't';
   'MISS_DISTANCE', 'f';
   'RELATIVE_SPEED', 'f';
   'RELATIVE_POSITION_R', 'f';
   'RELATIVE_POSITION_T', 'f';
   'RELATIVE_POSITION_N', 'f';
   'RELATIVE_VELOCITY_R', 'f';
   'RELATIVE_VELOCITY_T', 'f';
   'RELATIVE_VELOCITY_N', 'f';
   'Screening_Option', 's'; % CDM special COMMENT field
   'START_SCREEN_PERIOD', 't';
   'STOP_SCREEN_PERIOD', 't';
   'SCREEN_VOLUME_FRAME', 's';
   'SCREEN_VOLUME_SHAPE', 's';
   'SCREEN_VOLUME_X', 'f';
   'SCREEN_VOLUME_Y', 'f';
   'SCREEN_VOLUME_Z', 'f';
   'SCREEN_ENTRY_TIME', 't';
   'SCREEN_EXIT_TIME', 't';
   'COLLISION_PROBABILITY', 'f';
   'COLLISION_PROBABILITY_METHOD', 's'};

% List of object keywords
% CCSDS OBJECT keyword must be the first entry
% Column 1 is the keyword.  
% Column 2 is the value format.  's' - string, 'f' -  numeric, 't' - time
objParms = {
   'OBJECT', 's';
   'OBJECT_DESIGNATOR', 's';
   'CATALOG_NAME', 's';
   'OBJECT_NAME', 's';
   'INTERNATIONAL_DESIGNATOR', 's';
   'OBJECT_TYPE', 's';
   'OPERATOR_CONTACT_POSITION', 's';
   'OPERATOR_ORGANIZATION', 's';
   'OPERATOR_PHONE', 's';
   'OPERATOR_EMAIL', 's';
   'EPHEMERIS_NAME', 's';
   'COVARIANCE_METHOD', 's';
   'MANEUVERABLE', 's';
   'REF_FRAME', 's';
   'GRAVITY_MODEL', 's';
   'ATMOSPHERIC_MODEL', 's';
   'N_BODY_PERTURBATIONS', 's';
   'SOLAR_RAD_PRESSURE', 's';
   'EARTH_TIDES', 's';
   'INTRACK_THRUST', 's';
   'Covariance_Scale_Factor', 'f'; % CDM special COMMENT field
   'Exclusion_Volume_Radius', 'f'; % CDM special COMMENT field
   'TIME_LASTOB_START', 't';
   'TIME_LASTOB_END', 't';
   'RECOMMENDED_OD_SPAN', 'f';
   'ACTUAL_OD_SPAN', 'f';
   'OBS_AVAILABLE', 'f';
   'OBS_USED', 'f';
   'TRACKS_AVAILABLE', 'f';
   'TRACKS_USED', 'f';
   'RESIDUALS_ACCEPTED', 'f';
   'WEIGHTED_RMS', 'f';
   'Apogee', 'f';      % CDM special COMMENT field
   'Perigee', 'f';     % CDM special COMMENT field
   'Inclination', 'f'; % CDM special COMMENT field
   'AREA_PC', 'f';
   'MASS', 'f';
   'CD_AREA_OVER_MASS', 'f';
   'CR_AREA_OVER_MASS', 'f';
   'THRUST_ACCELERATION', 'f';
   'SEDR', 'f';
   'X', 'f';
   'Y', 'f';
   'Z', 'f';
   'X_DOT', 'f';
   'Y_DOT', 'f';
   'Z_DOT', 'f';
   'DCP_DENSITY_UNCERTAINTY', 'f'; % CDM special COMMENT field
   'DCP_SENSITIVITY_RTN_POS', 's'; % CDM special COMMENT field
   'DCP_SENSITIVITY_RTN_VEL', 's'; % CDM special COMMENT field
   'CR_R', 'f';
   'CT_R', 'f';
   'CT_T', 'f';
   'CN_R', 'f';
   'CN_T', 'f';
   'CN_N', 'f';
   'CRDOT_R', 'f';
   'CRDOT_T', 'f';
   'CRDOT_N', 'f';
   'CRDOT_RDOT', 'f';
   'CTDOT_R', 'f';
   'CTDOT_T', 'f';
   'CTDOT_N', 'f';
   'CTDOT_RDOT', 'f';
   'CTDOT_TDOT', 'f';
   'CNDOT_R', 'f';
   'CNDOT_T', 'f';
   'CNDOT_N', 'f';
   'CNDOT_RDOT', 'f';
   'CNDOT_TDOT', 'f';
   'CNDOT_NDOT', 'f';
   'CDRG_R', 'f';
   'CDRG_T', 'f';
   'CDRG_N', 'f';
   'CDRG_RDOT', 'f';
   'CDRG_TDOT', 'f';
   'CDRG_NDOT', 'f';
   'CDRG_DRG', 'f';
   'CSRP_R', 'f';
   'CSRP_T', 'f';
   'CSRP_N', 'f';
   'CSRP_RDOT', 'f';
   'CSRP_TDOT', 'f';
   'CSRP_NDOT', 'f';
   'CSRP_DRG', 'f';
   'CSRP_SRP', 'f'};

%% Open CDM file

[fid, errmsg] = fopen(filename,'rt');
if fid < 0
   fprintf('*** Error opening CDM file\n');
   fprintf('    %s\n', errmsg);
   fprintf('    File name: %s\n', strrep(filename,'\','\\'));
   Status = 1;
   return;
end

% Read all lines in the CDM file into data array
data = textscan(fid,'%s', 'Delimiter','\n');

% Close File
fclose(fid);

%% Process CDM text data array

% Extract the contents of the first cell to get the lines fo text
% and trim leading and trailing blanks from the lines of text
data = strtrim(data{1});

% Remove blank lines as represented by empty cells in the data array
idxValid = ~isempty_cell(data);
data = data(idxValid);

% Change special comment fields into effective keyword fields
idxComment = find(~isempty_cell(regexpi(data,'^COMMENT')));
NidxComment = numel(idxComment);
for ii = 1:NidxComment
    i = idxComment(ii);
    % disp(['Before: ' data{i}]);
    data{i} = strrep(data{i},'COMMENT Apogee Altitude','Apogee ');
    data{i} = strrep(data{i},'COMMENT Perigee Altitude','Perigee ');
    data{i} = strrep(data{i},'COMMENT Inclination','Inclination ');
    data{i} = strrep(data{i}, ...
        'COMMENT Screening Option', ...
        'Screening_Option ');
    data{i} = strrep(data{i}, ...
        'COMMENT Covariance Scale Factor', ...
        'Covariance_Scale_Factor ');
    data{i} = strrep(data{i}, ...
        'COMMENT Exclusion Volume Radius', ...
        'Exclusion_Volume_Radius ');
    data{i} = strrep(data{i}, ...
        'COMMENT DCP Density Forecast Uncertainty', ...
        'DCP_DENSITY_UNCERTAINTY ');
    data{i} = strrep(data{i}, ...
        'COMMENT DCP Sensitivity Vector RTN Pos', ...
        'DCP_SENSITIVITY_RTN_POS ');
    data{i} = strrep(data{i}, ...
        'COMMENT DCP Sensitivity Vector RTN Vel', ...
        'DCP_SENSITIVITY_RTN_VEL ');
    % disp(['After:  ' data{i}]);
end

% Remove comment lines (lines that start with the keword COMMENT)
idxValid = isempty_cell(regexpi(data,'^COMMENT'));
data = data(idxValid);

%% Get CDM version number, which must be in the first line
[keyword, value] = get_keyword_value(data{1});
% Verify thet the resulting keword is the CCSDS version keyword
word = headParms{1,1};
if ~strcmp(keyword,word)
    fprintf('*** Error reading CDM file.  First line does not contain %s\n', word);
    Status = 1;
   return
end
% Save version string to CDM header structure
cdmhead.(word) = value;

%% Find start of object 1 and object 2 data

% Finding the lines contain the OBJECT keyword
word = objParms{1,1};
idxObj = find(~isempty_cell(regexpi(data,['^' word '\s*='])));
if numel(idxObj) ~= 2
   fprintf('*** Error reading CDM file.  Expecting 2 objects.  Found %d objects\n',numel(idxObj));
    Status = 1;
   return
end
[keyword, value] = get_keyword_value(data{idxObj(1)});
if ~strcmpi(value,'OBJECT1')
    fprintf('*** Error reading CDM file.  First OBJECT (%s) not OBJECT1\n',value);
    Status = 1;
    return
end
[keyword, value] = get_keyword_value(data{idxObj(2)});
if ~strcmpi(value,'OBJECT2')
    fprintf('*** Error reading CDM file.  Second OBJECT (%s) not OBJECT2\n',value);
    Status = 1;
    return
end

%% Read remaining header and relative meta data parameters

% Last line of header and meta data is the line before the line containing
% the first OBJECT keyword
idxEnd = idxObj(1) - 1;

% Loop through each of the header/meta data keywords
for i = 2:size(headParms,1)
   word = headParms{i,1};  % Keyword
   fmt = headParms{i,2};   % Value format
   % Find the lines containing the current keyword at the beginning of the
   % line followed by a blank
   idx = find(~isempty_cell(regexpi(data(1:idxEnd),['^' word ' '])));
   
   % If more than one instance found, print a warning and use the last
   % instance
   if numel(idx) > 1
      fprintf('*** Warning -- More than one entry for keyword %s found\n',word);
      fprintf('    Using the last entry\n');
      idx = idx(end);
   end
   
   % If the keyword was found get and save its value
   if numel(idx) == 1
      [keyword, value] = get_keyword_value(data{idx}, fmt);
      if isempty(keyword) || isempty(value)
         fprintf('*** Warning -- Parsing of line for header/meta data keyword "%s" failed\n', word);
         fprintf('    Text: %s\n', data{idx});
      else
         cdmhead.(keyword) = value;
      end
   end
end % for each header/meta data keyword

%% Read OBJECT parameters for both objects
for iobj = 1:2
   % First data record for OBJECT i is the line after the ith occurance of
   % the OBJECT keyword
   idxStart = idxObj(iobj) + 1;
   
   % Last data record for OBJECT1 is the line before the line containing
   % the second OBJECT keyword
   % The data record for OBJECT2 is the last line of text
   if iobj == 1
      idxEnd = idxObj(iobj+1) - 1;
   else
      idxEnd = numel(data);
   end
   
   % Loop through each of the object data keywords
   for i = 2:size(objParms,1)
      word = objParms{i,1};   % Keyword
      fmt = objParms{i,2};    % Value format
      
      % Find the lines containing the current keyword within the data
      % reacord for the current object
      idx = find(~isempty_cell(regexpi(data(idxStart:idxEnd),['^' word ' '])));
      
      % If more than one instance found, print a warning and use the last
      % instance
      if numel(idx) > 1
         fprintf('*** Warning -- More than one entry for keyword %s found for OBJECT%d\n', word, iobj);
         fprintf('    Using the last entry\n');
         idx = idx(end);
      end
      
      % If the keyword was found get and save its value
      if numel(idx) == 1
         % disp(data{idxStart+idx-1})
         [keyword, value] = get_keyword_value(data{idxStart+idx-1}, fmt);
         % disp(keyword);
         % disp(value);
         if isempty(keyword) || isempty(value)
            fprintf('*** Warning -- Parsing of line for keyword "%s" failed for OBJECT%d\n', word, iobj);
            fprintf('    Text: %s\n', data{idxStart+idx-1});
         else
            cdmobj(iobj).(keyword) = value;
         end
      end
   end % for each object keyword
   
end % for each object

end  %Main function

function [keyword, value] = get_keyword_value(txt, fmt)

% Initialize output
keyword = '';
value = '';

if nargin < 2 || isempty(fmt)
   fmt = 's';
end

% Split "keyword = value [units]" text on equals sign
tokens = strsplit(txt,'=');
if numel(tokens) ~= 2
   return
end

% Keyword is first token with leading and trailing blanks removed
keyword = strtrim(tokens{1});

% Split second token ("value [units]") on "[" to remove unit labels if it
% exists
tokens = strsplit(tokens{2},'[');

% Value is first token with leading and trailing blanks removed
value = strtrim(tokens{1});

% COnvert to numeric format if specified
if strcmpi(fmt,'f')
   value = str2double(value);
elseif strcmpi(fmt,'t')
   timeFormat = getCcsdsTimeFormat(value);
   idx = strfind(timeFormat,'DDD');
   if ~isempty(idx)
      [DateNum,DateVec] = DOY2Date(str2num(value(idx(1):idx(1)+2)),str2num(value(1:4)));
      value = [num2str(DateVec(1),'%04d') '-' num2str(DateVec(2),'%02d') '-' num2str(DateVec(3),'%02d') 'T' value(idx(1)+4:end)];
   end
end

end % function get_keyword_value

function index_array = isempty_cell(input_cells)
    if ~iscell(input_cells)
       index_array = [];
       return
    end
    index_array = cellfun(@isempty,input_cells);
end

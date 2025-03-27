function [cdmhead, cdmobj, status] = read_cdm(filename, ignoreExtraFields, ignoreMissingEquals)
% Reads CCSDS CDM file
%
% Purpose:  Reads a CCSDS Conjunction Data Message (CDM) returning two
% data structures and a parse status. The first structure contains the CDM
% header and relative metadata information. The second structure contains
% the object 1 and 2 data.
%
% Invocation Method: [cdmhead, cdmobj] = read_cdm(filename)
%
% Argument            I/O  Description
% --------            ---  --------------------------------------------------
% filename             i   Fully qualified name of CDM file to read
% ignoreExtraFields    i   (Optional) When set to true, extra fields found
%                          in the CDM file will not cause the script to
%                          exit with errors.
%                          Default value = true
% ignoreMissingEquals  i   (Optional) When set to true, lines with no
%                          equals sign in them are ignored.
%                          Default value = true
% cdmhead              o   Structure containing the CDM header and relative
%                          metadata
% cdmobj(2)            o   Structure containing the object data for each of
%                          the two objects
% status               o   Variable indicating the status of parsing the
%                          CDM file: 0 = success, 1 = failure
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
% D. Hall       10/05/2020  Modified to handle odd characters within the
%                           DCP fields provided by CARA for SMAP event
% L. Baars      02/27/2023  Fixed relative pathing issue in addpath calls.
% L. Baars      04/21/2023  Added the ability to identify extra fields that
%                           appear within a CDM which aren't parsed by the
%                           reader. Also, allow for skipping of extra data
%                           at the beginning of a CDM before CCSDS_CDM_VERS
%                           is defined. Added the capability to move some
%                           fields from header to objects and vice-versa in
%                           order to fix some poorly formatted CDMs.
% L. Baars      02/14/2024  Added support for EFFECTIVE_HBR comment field.
% S. Es haghi   01/14/2025  Modified to read CDMs with a hidden
%                           initial character in their CCSDS keyword token

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'../TimeTransformations')); addpath(s.path);
    pathsAdded = true;
end

if nargin == 1
    ignoreExtraFields = true;
    ignoreMissingEquals = true;
elseif nargin == 2
    ignoreMissingEquals = true;
elseif nargin ~= 3
    error('Invalid number of arguments passed in!');
end

% Initialize output
cdmhead = [];
cdmobj  = [];
status = 0;

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
   'COLLISION_PROBABILITY_METHOD', 's';
   'Screened_With', 's';               % CDM special COMMENT field, object data that sometimes appears in the header
   'SCREENING_DATA_SOURCE', 's';       % CDM special COMMENT field, object data that sometimes appears in the header
   'HBR', 'f'};                        % CDM special COMMENT field, header data that sometimes appears in the object

headerMap = containers.Map('KeyType','char','ValueType','logical');
for i = 1:size(headParms,1)
    headerMap(headParms{i,1}) = true;
end

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
   'OD_DATA_SOURCE', 's';
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
   'CSRP_SRP', 'f';
   'Screened_With', 's';               % CDM special COMMENT field, object data that sometimes appears in the header
   'SCREENING_DATA_SOURCE', 's';       % CDM special COMMENT field, object data that sometimes appears in the header
   'HBR', 'f'};                        % CDM special COMMENT field, header data that sometimes appears in the object

objectMap = containers.Map('KeyType','char','ValueType','logical');
for i = 1:size(objParms,1)
    objectMap(objParms{i,1}) = true;
end

%% Open CDM file

% Open CDM file
[fid, errmsg] = fopen(filename,'rt');
% =========================================================================
% Open file using UTF-8 native format to read CARA-produced CDMs
% [fid, errmsg] = fopen(filename,'r','n','UTF-8');
% =========================================================================

if fid < 0
   fprintf('*** Error opening CDM file\n');
   fprintf('    %s\n', errmsg);
   fprintf('    File name: %s\n', strrep(filename,'\','\\'));
   status = 1;
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

% Remove lines with no equals sign, but only if ignoreMissingEquals is set
% to true
if ignoreMissingEquals
    keepIdx = ~isempty_cell(regexpi(data,'='));
    data = data(keepIdx);
else
    errorIdx = isempty_cell(regexpi(data,'='));
    if sum(errorIdx) > 0
        fprintf('*** Error reading CDM file. Found %d lines without equals signs\n', sum(errorIdx));
        fprintf('    Text:\n');
        data = data(errorIdx);
        for i = 1:length(data)
            fprintf('      %s\n',data{i});
        end
        status = 1;
        return;
    end
end

% Change special comment fields into effective keyword fields
idxComment = find(~isempty_cell(regexpi(data,'^COMMENT')));
NidxComment = numel(idxComment);
for ii = 1:NidxComment
    i = idxComment(ii);
    % disp(['Before: ' data{i}]);
    data{i} = strrep(data{i},'COMMENT Apogee Altitude','Apogee ');
    data{i} = strrep(data{i},'COMMENT Perigee Altitude','Perigee ');
    data{i} = strrep(data{i},'COMMENT Inclination','Inclination ');
    % Process Screening_Option comments (JSPOC and CARA formats)
    data{i} = strrep(data{i}, ...
        'COMMENT Screening Option', ...
        'Screening_Option ');
    data{i} = strrep(data{i}, ...
        'COMMENT SCREENING_OPTION', ...
        'Screening_Option ');
    % Process Covariance_Scale_Factor comments (JSPOC and CARA formats)
    data{i} = strrep(data{i}, ...
        'COMMENT Covariance Scale Factor', ...
        'Covariance_Scale_Factor ');
    data{i} = strrep(data{i}, ...
        'COMMENT COVARIANCE_SCALE_FACTOR', ...
        'Covariance_Scale_Factor ');
    % Process Exclusion_Volume_Radius comments (JSPOC and CARA formats)
    data{i} = strrep(data{i}, ...
        'COMMENT Exclusion Volume Radius', ...
        'Exclusion_Volume_Radius ');
    data{i} = strrep(data{i}, ...
        'COMMENT EXCLUSION_VOLUME_RADIUS', ...
        'Exclusion_Volume_Radius ');
    % Substitute odd characters found in CARA-produced CDMs for a single
    % space character. Needed for SMAP Ops Request CDMs 10/5/2020.
    if ~isempty(strfind(data{i},'DCP'))
        if ~isempty(strfind(data{i},char(160))) || ...
           ~isempty(strfind(data{i},char(194)))
            disp('WARNING: Replacing extraneous character(s) in DCP line');
            disp(['Before replacement: ' data{i}])
            data{i} = strrep(data{i},char(194),'');
            data{i} = strrep(data{i},char(160),' ');
            disp(['After  replacement: ' data{i}])
        end
    end
    % OD Data source
    data{i} = strrep(data{i}, ...
        'COMMENT OD_DATA_SOURCE', ...
        'OD_DATA_SOURCE');
    % Covariance cross-correction params
    data{i} = strrep(data{i}, ...
        'COMMENT DCP Density Forecast Uncertainty', ...
        'DCP_DENSITY_UNCERTAINTY ');
    data{i} = strrep(data{i}, ...
        'COMMENT DCP Sensitivity Vector RTN Pos', ...
        'DCP_SENSITIVITY_RTN_POS ');
    data{i} = strrep(data{i}, ...
        'COMMENT DCP Sensitivity Vector RTN Vel', ...
        'DCP_SENSITIVITY_RTN_VEL ');
    data{i} = strrep(data{i}, ...
        'COMMENT SCREENING_DATA_SOURCE', ...
        'SCREENING_DATA_SOURCE ');
    % Hard-body radius
    data{i} = strrep(data{i}, ...
        'COMMENT Operator Hard Body Radius', ...
        'HBR ');
    data{i} = strrep(data{i}, ...
        'COMMENT HBR', ...
        'HBR ');
    data{i} = strrep(data{i}, ...
        'COMMENT EFFECTIVE_HBR', ...
        'HBR ');
    % Screened with special processing, with and without the '=' sign
    screenedWithComment = 'COMMENT Screened with';
    if ~isempty(regexpi(data{i},screenedWithComment)) && ...
            isempty(regexpi(data{i},[screenedWithComment ' *=']))
        data{i} = strrep(data{i}, ...
            screenedWithComment, ...
            'Screened_With =');
    else
        data{i} = strrep(data{i}, ...
            screenedWithComment, ...
            'Screened_With');
    end
    % disp(['After:  ' data{i}]);
end

% Find comment lines that haven't been defined and could have a value
idxUndefComments = ~isempty_cell(regexpi(data,'^COMMENT.*=')) & isempty_cell(regexpi(data,'^COMMENT *='));
if sum(idxUndefComments) > 0
    tempData = data(idxUndefComments);
    if ignoreExtraFields
        fprintf('*** Warning -- Found %d undefined comment fields:\n', sum(idxUndefComments));
    else
        fprintf('*** Error reading CDM file. Found %d undefined comment fields:\n', sum(idxUndefComments));
    end
    for i = 1:numel(tempData)
        fprintf('    Text: %s\n',tempData{i});
    end
    if ~ignoreExtraFields
        status = 1;
        return;
    end
end

% Remove remaing comment lines (lines that start with the keword COMMENT
% and didn't have an equals sign identifying a name/value pair)
idxValid = isempty_cell(regexpi(data,'^COMMENT'));
data = data(idxValid);

%% Get CDM version number, which will start the data section of the CDM
% Verify thet the resulting keword is the CCSDS version keyword
word = headParms{1,1};
% Loop until we find the CCSDS keyword or we get to the end of the file
keyword = '';
headStartIdx = 0;
while ~strcmp(keyword,word) && headStartIdx < numel(data)
    headStartIdx = headStartIdx + 1;
    [keyword, value] = get_keyword_value(data{headStartIdx});
    if length(keyword) == 15; keyword = keyword(2:end); end % Deletes the hidden character in special case CDMs
end
% Check if the CCSDS keyword was found
if ~strcmp(keyword,word)
    fprintf('*** Error reading CDM file. The %s token was not found\n', word);
    status = 1;
    return
end
if headStartIdx ~= 1
    fprintf('*** Warning -- Ignored %d lines before the %s token was found\n', headStartIdx-1, word);
end
% Save version string to CDM header structure
cdmhead.(word) = value;

%% Find start of object 1 and object 2 data

% Finding the lines contain the OBJECT keyword
word = objParms{1,1};
idxObj = find(~isempty_cell(regexpi(data,['^' word '\s*='])));
if numel(idxObj) ~= 2
   fprintf('*** Error reading CDM file. Expecting 2 objects. Found %d objects\n',numel(idxObj));
    status = 1;
   return
end
[~, value] = get_keyword_value(data{idxObj(1)});
if ~strcmpi(value,'OBJECT1')
    fprintf('*** Error reading CDM file. First OBJECT (%s) not OBJECT1\n',value);
    status = 1;
    return
end
[~, value] = get_keyword_value(data{idxObj(2)});
if ~strcmpi(value,'OBJECT2')
    fprintf('*** Error reading CDM file. Second OBJECT (%s) not OBJECT2\n',value);
    status = 1;
    return
end

%% Read remaining header and relative meta data parameters

% Last line of header and meta data is the line before the line containing
% the first OBJECT keyword
idxEnd = idxObj(1) - 1;

% Loop through each of the header/meta data keywords
tempData = data(headStartIdx+1:idxEnd);
for i = 2:size(headParms,1)
   word = headParms{i,1};  % Keyword
   fmt = headParms{i,2};   % Value format
   % Find the lines containing the current keyword at the beginning of the
   % line followed by a blank
   idx = find(~isempty_cell(regexpi(tempData,['^' word ' '])));
   
   % If more than one instance found, print a warning and use the last
   % instance
   if numel(idx) > 1
      fprintf('*** Warning -- More than one entry for keyword %s found\n',word);
      fprintf('    Using the last entry\n');
      idx = idx(end);
   end
   
   % If the keyword was found get and save its value
   if numel(idx) == 1
      [keyword, value] = get_keyword_value(tempData{idx}, fmt);
      if isempty(keyword) || isempty(value)
         fprintf('*** Warning -- Parsing of line for header/meta data keyword "%s" failed\n', word);
         fprintf('    Text: %s\n', tempData{idx});
      else
         cdmhead.(keyword) = value;
      end
   end
end % for each header/meta data keyword

%% Check for extra keywords in the header which aren't in our definition
numMissing = 0;
for i = headStartIdx+1:idxEnd
    keyword = get_keyword_value(data{i});
    if ~isKey(headerMap,keyword)
        numMissing = numMissing + 1;
        if numMissing == 1
            if ignoreExtraFields
                fprintf('*** Warning -- At least one header keyword was not found in the header map:\n');
            else
                fprintf('*** Error reading CDM file. At least one header keyword was not found in the header map:\n');
            end
        end
        fprintf('    Text: %s\n',data{i});
    end
end
if numMissing > 0 && ~ignoreExtraFields
    status = 1;
    return;
end

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
   
   % Get the keywords found between idxStart and idxEnd
   dataMap = containers.Map('KeyType','char','ValueType','double');
   numMissing = 0;
   for i = idxStart:idxEnd
       dataline = data{i};
       if isempty(dataline)
           continue;
       end
       firstSpaceIdx = find(isspace(dataline),1);
       if isempty(firstSpaceIdx)
           stopIdx = length(data{i});
       elseif firstSpaceIdx > 1
           stopIdx = firstSpaceIdx - 1;
       else
           continue;
       end
       keyword = dataline(1:stopIdx);
       % Check for extra keywords in the object which aren't in our definition
       if ~isKey(objectMap,keyword)
           numMissing = numMissing + 1;
            if numMissing == 1
                if ignoreExtraFields
                    fprintf('*** Warning -- At least one object keyword was not found in the object map:\n');
                else
                    fprintf('*** Error reading CDM file. At least one object keyword was not found in the object map:\n');
                end
            end
            fprintf('    Text: %s\n',dataline);
       end
       % Add the index of the current keyword to our dataMap
       if ~isKey(dataMap,keyword)
           dataMap(keyword) = i;
       else
           % Use negative numbers to indicate that more than one instance
           % of a keyword has been found
           dataMap(keyword) = -i;
       end
   end
   if numMissing > 0 && ~ignoreExtraFields
        status = 1;
        return;
   end
   
   % Loop through each of the object data keywords
   for i = 2:size(objParms,1)
      word = objParms{i,1};   % Keyword
      fmt = objParms{i,2};    % Value format
      
      % Find the lines containing the current keyword within the data
      % record for the current object
      if isKey(dataMap,word)
          idx = dataMap(word);
          if idx < 0
             % If more than one instance found, print a warning and use the
             % last instance
             fprintf('*** Warning -- More than one entry for keyword %s found for OBJECT%d\n', word, iobj);
             fprintf('    Using the last entry\n');
             idx = -idx;
          end
      else
          idx = [];
      end
      
      % If the keyword was found get and save its value
      if numel(idx) == 1
         % disp(data{idx})
         [keyword, value] = get_keyword_value(data{idx}, fmt);
         % disp(keyword);
         % disp(value);
         if isempty(keyword) || isempty(value)
            fprintf('*** Warning -- Parsing of line for keyword "%s" failed for OBJECT%d\n', word, iobj);
            fprintf('    Text: %s\n', data{idx});
         else
            cdmobj(iobj).(keyword) = value; %#ok<AGROW>
         end
      end
   end % for each object keyword
   
end % for each object

%% Cleanup object values that sometimes appear in the header
headerToObj = {'Screened_With','SCREENING_DATA_SOURCE'};
for i = 1:length(headerToObj)
    keyword = headerToObj{i};
    if isfield(cdmhead,keyword) && isfield(cdmobj,keyword)
        cdmobj(2).(keyword) = cdmobj(1).(keyword);
        cdmobj(1).(keyword) = cdmhead.(keyword);
        cdmhead = rmfield(cdmhead,keyword);
    end
end

%% Cleanup header values that sometimes appear in the object
objToHeader = {'HBR'};
for i = 1:length(objToHeader)
    keyword = objToHeader{i};
    if isfield(cdmobj,keyword) && ~isfield(cdmhead,keyword)
        cdmhead.(keyword) = cdmobj(1).(keyword);
        cdmobj = rmfield(cdmobj,keyword);
    end
end

end  %Main function

function [keyword, value] = get_keyword_value(txt, fmt)

% Initialize output
keyword = '';
value = '';

if nargin < 2 || isempty(fmt)
   fmt = 's';
end

% Split "keyword = value [units]" text on equals sign
idxs = strfind(txt,'=');
if length(idxs) ~= 1
    return;
end
tokens = cell(2,1);
tokens{1} = txt(1:idxs-1);
tokens{2} = txt(idxs+1:end);

% Keyword is first token with leading and trailing blanks removed
keyword = strtrim(tokens{1});

% Split second token ("value [units]") on "[" to remove unit labels if it
% exists
value = tokens{2};
idxs = strfind(value,'[');
if length(idxs) >= 1
    value = value(1:idxs(1)-1);
end
value = strtrim(value);

% Convert to numeric format if specified
if strcmpi(fmt,'f')
   value = str2double(value);
elseif strcmpi(fmt,'t')
   timeFormat = getCcsdsTimeFormat(value);
   if isempty(timeFormat)
       value = '';
   else
       idx = strfind(timeFormat,'DDD');
       if ~isempty(idx)
          [~,DateVec] = DOY2Date(str2double(value(idx(1):idx(1)+2)),str2double(value(1:4)));
          value = [num2str(DateVec(1),'%04d') '-' num2str(DateVec(2),'%02d') '-' num2str(DateVec(3),'%02d') 'T' value(idx(1)+4:end)];
       end
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

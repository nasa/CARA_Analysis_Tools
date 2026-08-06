function [CCSDSData,Status,data_text] = CCSDSParser(CCSDSEphemFile)
% CCSDSParser - Parse state and covariance information from a CCSDS
%               ephemeris file.
%
% Syntax: [CCSDSData, Status, data_text] = CCSDSParser(CCSDSEphemFile);
%
% =========================================================================
%
% Copyright (c) 2015-2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    CCSDSEphemFile  -  CCSDS file location and name
%
% =========================================================================
%
% Output:
%
%    CCSDSData       -  Structure with fields:
%                        .Epoch - Matlab datenum             [Nx1]
%                        .State - Position and velocity      [Nx6]
%                        .Cov   - Covariance                 [Nx36]
%
%    Status          -  Success flag (0 = passed, 1 = failed,
%                       2 = no covariance data)
%
%    data_text       -  File data from textscan
%                       (without empty lines)
%
% =========================================================================
%
% Initial version: Feb 2015;  Latest update: Jul 2026
%
% ----------------- BEGIN CODE -----------------

% Initialize output parameters
CCSDSData = struct('Epoch',[],'State',[],'Cov',[]);
Status = 0;

try
    %% Read CCSDS ephemeris/cov file
    [FID, msg]  = fopen(CCSDSEphemFile,'rt');
    if ~isempty(msg)
        fprintf('Error opening CCSDS Ephemeris file: %s\n     %s\n',CCSDSEphemFile, msg);
        Status = 1;
        return
    end
    data = textscan(FID,'%s','Delimiter','');
    fclose(FID);
    data = strtrim(data{1});
    
    % Remove blank lines
    data = data(~isempty_cell(data));
    
    % Save original text data
    data_text = data;
    
    % Eliminate comment lines
    notComment = isempty_cell(regexpi(data,'COMMENT'));    
    data = data(notComment);
    if any(~notComment)
        warning(['Eliminating ' num2str(sum(~notComment)) ' COMMENT lines before data parsing']);
    end
    
    % Get size of CCSDS file
    rowNum = numel(data);
    
    timeIdx = find(~isempty_cell(regexpi(data,'TIME_SYSTEM')));
    timeSys = textscan(char(data(timeIdx)),'%*s %*s %s');
    
    % Determine indices for start and stop of ephemeris points
    startEphemIdx = find(~isempty_cell(regexpi(data,'META_STOP')));
    startEphemIdx = startEphemIdx(1) + 1;
    stopEphemIdx = find(~isempty_cell(regexpi(data,'COVARIANCE_START')));
    if ~isempty(stopEphemIdx)
        stopEphemIdx = stopEphemIdx(1) - 1;
    else
        stopEphemIdx = rowNum;
    end
    
    %% Parse state data
    npts = stopEphemIdx - startEphemIdx + 1;
    % Must append blanks to data to ensure textscan works
    blanks = repmat(' ', npts, 1);
    d = [char(data(startEphemIdx:stopEphemIdx)), blanks];
    
    % Check first line to get number of columns of state data
    d1 = textscan(d(1,:)','%s %f %f %f %f %f %f %f %f %f');
    a = [d1{2:end}];
    ncol = numel(a);
    % Set read format based on number of columns of state data
    fmt = ['%s' repmat(' %f',1,ncol)];
    
    % Parse all lines of state data
    stateData = textscan(d',fmt);
    
    % Get epoch time format to speed up conversion
    % Remove T from epoch time for Matlab
    Epoch = strrep(stateData{1},'T',' ');
    % Get number of digits in fractions of seconds in epoch time
    tstr = Epoch{1};
    idx = strfind(tstr,'.');
    if isempty(idx)
        fstr = '';
    else
        nfrac = numel(tstr) - idx;
        nfrac = min(nfrac,6); % Up to microsecond if applicable
        fstr = ['.' repmat('S',1,nfrac)];
    end
    time_formatMain = ['uuuu-MM-dd HH:mm:ss' fstr];
    time_formatBackup = ['uuuu-DDD HH:mm:ss' fstr];
    
    % Convert Epoch time to Matlab date number and save to CCSDSData structure
    try
        CCSDSData.Epoch = datenum(datetime(Epoch ,'InputFormat', time_formatMain));
    catch
        CCSDSData.Epoch = datenum(datetime(Epoch ,'InputFormat', time_formatBackup));
    end
    
    % Convert time to UTC with proper amount of leapseconds
    if strcmpi(timeSys{1},'TAI')
        Status = 1;
        fprintf('CCSDS Ephemeris file is in TAI convert to UTC: %s\n',CCSDSEphemFile);
        return
    end
    
    % Concantenate first 6 columns of state data and save to CCSDSData
    % structure
    State = [stateData{2}, stateData{3}, stateData{4}, stateData{5},...
        stateData{6}, stateData{7}];
    CCSDSData.State = State;
    
    % If no covariance data, return with a covariance matrix of zeros
    if (stopEphemIdx == rowNum)
        Status = 2;
        covdata = zeros(npts,36); %CCSDSParser add upper triangle to cov before returning, so just doing a 6x6
        CCSDSData(1,1).Cov   = covdata;
        return
    end
    
    %% Parse covariance data
    % Initialize format for each line of covariance data
    covfmt = {'%f'; '%f %f'; '%f %f %f'; '%f %f %f %f'; '%f %f %f %f %f'; '%f %f %f %f %f %f'};
    
    % Get locations of lines containing covariance epoch time and
    % covariance reference frame parameters
    idxepoch = find(~isempty_cell(regexpi(data,'EPOCH')));
    idxframe = find(~isempty_cell(regexpi(data,'COV_REF_FRAME')));
    % Set locations of the first line of each set of covariance data
    if isempty(idxframe)
        % No reference frame parameter, 
        % first line of covariance data is after the epoch data
        idxcov1 = idxepoch + 1;
    else
        % Reference frame parameter defined,
        % first line of covariance data is after the reference frame data
        if all(idxframe>idxepoch)
            idxcov1 = idxframe + 1;
        else
            idxcov1 = idxepoch + 1;
        end
    end

    % Get number of lines of covariance elements
    if numel(idxepoch) > 1
        if isempty(idxframe)
        ncov = idxepoch(2) - idxcov1(1);
        else
            if all(idxframe>idxepoch)
                ncov = idxepoch(2) - idxcov1(1);
            else
                ncov = idxframe(2) - idxcov1(1);
            end
        end
    else
        ncov = rowNum  - idxcov1(1);
    end
    % Get number of covariance sets
    npts = numel(idxcov1);
    % Initialize blanks to add to end of each line of data for proper textscan
    blanks = repmat(' ', npts, 1);
    % Initialize covariance data array
    covdata = zeros(npts,ncov*(ncov+1)/2);
    
    % Read the first line of covariances for each set followed by the second
    % line, then the third line, etc. and append the data to the proper
    % columns in the covariance data array
    idx = 1;
    for i = 1:ncov
        txt = [char(data(idxcov1+i-1)) blanks];
        d = textscan(txt',covfmt{i});
        for j = 1:i
            covdata(:,idx+j-1) = d{j};
        end
        idx = idx + i;
    end
    % Convert the lower trianglular covariance matrix to the full covariance
    % matrix
    currCov  = triu(ones(ncov));
    currCov(currCov~=0) = 1:ncov*(ncov+1)/2;
    currCov  = reshape(currCov + triu(currCov,1)', 1, ncov*ncov);
    covdata = covdata(:,currCov);
    
    % Save the covariances to the CCSDSData structure
    CCSDSData(1,1).Cov   = covdata;
    
catch Me
    Status = 1;
    fprintf('Error in CCSDS Parser for ephemeris file: %s\n     %s\n',CCSDSEphemFile,Me.message);
end

return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer     |     Date    | Description
% ---------------------------------------------------
% D. Plakalovic | 2015-Feb-09 | Created
% R. Coon       | 2016-Aug-29 | Sped up processing, added checks for number
%               |             | of covariance elements, epoch time format, 
%               |             | and presence or absence of covariance 
%               |             | reference frame parameter
% L. Johnson    | 2016-Dec-05 | Modified error handling. Populate blank
%               |             | covariance with zeros.
% J. Halpin     | 2026-Apr-22 | Reformatted header to match code standards.
%               |             | Added a footer.
% S. Es haghi   | 2026-Jun-25 | Add capability to handle epochs down to
%               |             | seconds
% S. Es haghi   | 2026-Jul-01 | Add capability to handle days of year
%               |             | epochs and cases where covariance reference
%               |             | frame and epochs lines are in the reverse
%               |             | order
% =========================================================================
%
% Copyright (c) 2015-2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
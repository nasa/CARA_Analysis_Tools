function [discosData] = GetDiscosData(discosFile, connParams, excludeBackspaceInLogOutputs)
% GetDiscosData - Get DISCOS Data from DISCOS website
%
% Syntax:
%
%   [discosData] = GetDiscosData(discosFile);
%   [discosData] = GetDiscosData(discosFile, connParams);
%   [discosData] = GetDiscosData(discosFile, connParams,excludeBackspaceInLogOutputs);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Inputs:
%
%   discosFile - Name of the file where downloaded DISCOS data should be
%                stored. If the file exists and has data in it, then the
%                data within the file is returned back to the calling
%                function. Otherwise, the script will download the data
%                from the DISCOS website, store that data into the
%                discosFile, and return the data back to the calling
%                function.
%
%   connParams - (Optional) Connection parameters associated with
%                connecting to the DISCOS website. If this parameter is not
%                passed in, then it is assumed that the connection
%                parameters are stored in the file connect_params.m.
%
%   excludeBackspaceInLogOutputs - (Optional) When set to true, the log
%                                  outputs produced by the tool will
%                                  produce a running log. The default for
%                                  this parameter is false, which will
%                                  prevent the log data from flooding the
%                                  screen.
%
% =========================================================================
%
% Outputs:
%
%   discosData - An nx14 table containing information about each object
%                retrieved from the DISCOS database. Field included are:
%     ObjectID - Satellite Catalog NORAD catalog ID
%     ObjectName - Common name for the object
%     ObjectClass - Type of object as classified by DISCOS (e.g. Rocket
%                   Body, Payload, etc.)
%     ObjectShape - Object shape as classified by DISCOS (e.g. Cyl, Sphere,
%                   Box, etc.)
%     Width - DISCOS defined width, in meters (value is populated depending
%             on the shape of the object)
%     Height - DISCOS defined height, in meters (value is populated
%              depending on the shape of the object)
%     Depth - DISCOS defined depth, in meters (value is populated depending
%             on the shape of the object)
%     Diameter - DISCOS defined diameter, in meters (value is populated
%                depending on the shape of the object)
%     Span - DISCOS defined span, in meters. This is the maximum extent of
%            any single dimension of the object.
%     HBR - CARA calculated Hard Body Radius (HBR) based on the DISCOS
%           provided dimension and shape parameters.
%     MinCrossSection - DISCOS defined minimum cross sectional area, in m^2
%     MaxCrossSection - DISCOS defined maximum cross sectional area, in m^2
%     AvgCrossSection - DISCOS defined average cross sectional area, in m^2
%     Mass - DISCOS defined mass, in kg
%
% =========================================================================
%
% Initial version: Sep 2021;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    % Checka nd set default arguments
    if nargin == 1
        connParams = connect_params;
        excludeBackspaceInLogOutputs = false;
    elseif nargin == 2
        excludeBackspaceInLogOutputs = false;
    elseif nargin ~= 3
        error('Invalid number of parameters passed in!');
    end
    
    % If the file exists, is not empty, and can be read then don't download
    % it again
    s = dir(discosFile);
    if ~isempty(s) && s.bytes > 0
        try
            discosData = readtable(discosFile);
            % If we got here, then we were able to read the file, no need
            % to download it again
            return;
        catch
        end
    end
    
    % Open and close the output file just to make sure we can write to it
    % before we query the DISCOS website. The download from the website
    % takes a long time to complete and we want to make sure that we can
    % actually write to the output file before we run any of these queries.
    try
        fid = fopen(discosFile,'w');
        fclose(fid);
    catch me
        error(['Could not open ' discosFile ' for writing: ' me.message]);
    end
    
    % Set URL and login credentials
    connidx = GetConnectionIdx('DISCOS', connParams);
    connHost = connParams.conn(connidx).conn_host;
    baseURL = ['https://' connHost];
    secKey = GetSecKey;
    discosToken = GetPassword('DISCOS',secKey,connParams);
    
    % Set query parameters
    pageSize = 100; % maxes out at 100
    page = 1;       % get first page
    queryString = [baseURL '/api/objects?filter=ne(satno,null)&sort=satno&page%%5Bsize%%5D=%d&page%%5Bnumber%%5D=%d'];
    
    % Set the options for passing the DISCOS token and the DISCOS API
    % version
    options = weboptions('HeaderFields',{'Authorization' ['Bearer ' discosToken];'DiscosWeb-Api-Version' '2'});
    
    % Begin the download process
    startText = 'Downloading DISCOS data...';
    disp([startText 'page 1 of ?']);
    URL = sprintf(queryString,pageSize,page);
    % The next line actually downloads the first page of data
    tempData = webread(URL,options);
    numDel = 12;
    if isstruct(tempData)
        % Get the total page number from the meta data, we can't get the
        % total page number otherwise
        totalPages = tempData.meta.pagination.totalPages;
        % Calculate the total number of elements we'll retrieve and setup
        % the memory to hold all of the data, we don't want these arrays
        % increasing in size on each iteration
        maxElems = totalPages * pageSize;
        ObjectID = nan(maxElems,1);
        ObjectName = cell(maxElems,1);
        ObjectClass = cell(maxElems,1);
        ObjectShape = cell(maxElems,1);
        Width = nan(maxElems,1);
        Height = nan(maxElems,1);
        Depth = nan(maxElems,1);
        Diameter = nan(maxElems,1);
        Span = nan(maxElems,1);
        MinCrossSection = nan(maxElems,1);
        MaxCrossSection = nan(maxElems,1);
        AvgCrossSection = nan(maxElems,1);
        Mass = nan(maxElems,1);
        idx = 1;
        % Continue with the download process until we've retrieved all
        % pages
        while page <= totalPages
            % Get the actual number of elements on the current page
            numInStruct = length(tempData.data);
            numInLoop = min(numInStruct,pageSize);
            % Loop through all of the elements on the current page
            for i = 1:numInLoop
                % Save off the values for each element on the page
                ObjectID(idx) = GetValidValue(tempData.data(i).attributes.satno);
                ObjectName{idx} = tempData.data(i).attributes.name;
                ObjectClass{idx} = tempData.data(i).attributes.objectClass;
                ObjectShape{idx} = tempData.data(i).attributes.shape;
                Width(idx) = GetValidValue(tempData.data(i).attributes.width);
                Height(idx) = GetValidValue(tempData.data(i).attributes.height);
                Depth(idx) = GetValidValue(tempData.data(i).attributes.depth);
                Diameter(idx) = GetValidValue(tempData.data(i).attributes.diameter);
                Span(idx) = GetValidValue(tempData.data(i).attributes.span);
                MinCrossSection(idx) = GetValidValue(tempData.data(i).attributes.xSectMin);
                MaxCrossSection(idx) = GetValidValue(tempData.data(i).attributes.xSectMax);
                AvgCrossSection(idx) = GetValidValue(tempData.data(i).attributes.xSectAvg);
                Mass(idx) = GetValidValue(tempData.data(i).attributes.mass);
                % Increase the index in our allocated memory so we don't
                % overwrite what we just saved
                idx = idx + 1;
            end
            % Increase the page number and see if we need to continue
            page = page + 1;
            if page <= totalPages
                % Some special handling for logging to an output file vs.
                % writing to screen without filling up the output window
                if excludeBackspaceInLogOutputs
                    bs = startText;
                else
                    bs = repmat('\b',1,numDel);
                end
                pageTxt = ['page ' num2str(page) ' of ' num2str(totalPages)];
                fprintf([bs '%s\n'], pageTxt);
                numDel = length(pageTxt) + 1;
                % Setup the URL for the next page query
                URL = sprintf(queryString,pageSize,page);
                try
                    % Try to download the next page
                    tempData = webread(URL,options);
                catch
                    % We'll get a download error when we've requested for
                    % too much data in too short of an amount of time. This
                    % code will pause for 60 seconds and then attempt to
                    % redownload the page.
                    dispTxt = ' (pausing to prevent "TOO MANY REQUESTS" response)';
                    if excludeBackspaceInLogOutputs
                        fprintf('%s\n',dispTxt);
                    else
                        fprintf('\b%s\n',dispTxt);
                    end
                    numDel = numDel + length(dispTxt);
                    pause(60);
                    tempData = webread(URL,options);
                end
                % If no data was found, then we had problems downloading
                % the data and just need to exit
                if ~isstruct(tempData)
                    error(['Could not read page ' num2str(page) ' of ' num2str(totalPages) ' from ' connHost]);
                end
            end
        end
    else
        error(['Could not read the first page from ' connHost]);
    end
    if excludeBackspaceInLogOutputs
        fprintf('success\n');
    else
        bs = repmat('\b',1,numDel);
        fprintf([bs 'success\n']);
    end
    
    % Calculate the HBR
    numElems = idx - 1;
    HBR = nan(maxElems, 1);
    for i = 1:numElems
        shape = lower(ObjectShape{i});
        if isempty(shape)
            continue;
        end
        x = Width(i);
        y = Height(i);
        z = Depth(i);
        diam = Diameter(i);
        span = Span(i);
        if ~startsWith(shape,'sphere')
            if ~isnan(x) && ~isnan(y) && ~isnan(z) && ~isnan(span)
                % We assume that anything that has all four fields filled
                % in will be enclosed either by the max extent defined by
                % the x, y, and z dimensions or the span. Estimate the HBR
                % from the max of these two values.
                val1 = sqrt(x^2 + y^2 + z^2) / 2;
                val2 = span / 2;
                HBR(i) = max(val1, val2);
                if strcmp(shape,'box') && val2 > val1
                    ObjectShape{i} = 'Box with extra span';
                end
            elseif (contains(shape,'cone') || contains(shape,'cyl')) && ~isnan(y) && ~isnan(diam)
                if contains(shape,'cone')
                    % The max extent of a cone is from the tip of the cone
                    % (at diam/2) to the edge of the circle at the bottom
                    % of the cone
                    val1 = sqrt(y^2 + (diam/2)^2) / 2;
                else
                    % The max extent of a cylinder is from the edge of the
                    % circle at the top of the cylinder to the opposite
                    % edge of the circle at the bottom of the cylinder
                    val1 = sqrt(y^2 + diam^2) / 2;
                end
                if ~isnan(span)
                    % Use the larger of the two values, the max extent
                    % calculated above, or the span
                    HBR(i) = max(val1, span / 2);
                else
                    HBR(i) = val1;
                end
            elseif ~isnan(x) && ~isnan(y) && ~isnan(z)
                % If span isn't defined and we haven't fallen into one of
                % the other cases, use x, y, and z to estimate the HBR
                HBR(i) = sqrt(x^2 + y^2 + z^2) / 2;
            elseif ~isnan(span)
                % If span is defined and we haven't calculated an HBR, use
                % the span
                HBR(i) = span / 2;
            end
        else
            % For spherical and ellipsoidal objects, the HBR of the
            % circumscribing sphere is half of the largest dimension
            HBR(i) = max([x y z diam span]) / 2;
        end
    end
    
    % Find indexes where we should cleanup the data in case DISCOS didn't
    % provide some data for us
    removeIdx = isnan(ObjectID) | (isnan(HBR) & isnan(Mass));
    % Find objects in reserved ranges, remove those as well
    removeIdx = removeIdx | (ObjectID >= 70000 & ObjectID <= 99999) | ...
        (ObjectID >= 700000000 & ObjectID <= 999999999);
    % Create a table based on the data that was read and calculated
    discosData = table(ObjectID,ObjectName,ObjectClass,ObjectShape,Width,Height,Depth,Diameter,Span,HBR,MinCrossSection,MaxCrossSection,AvgCrossSection,Mass);
    % Cleanup the data
    discosData(removeIdx,:) = [];
    writetable(discosData,discosFile);
    disp(['Saved DISCOS data into file: ' discosFile]);
end

% Sets empty or zero values into NaN values. A size of zero is not valid in
% the context of sizing information, so we override it with a NaN.
function [validValue] = GetValidValue(value)
    if isempty(value) || value == 0
        validValue = nan;
    else
        validValue = value;
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2021-Sep-08 | Initial Development.
% D. Hall        | 2022-May-23 | Expanded output table to include other
%                |             | fields in addition to the HBR and Mass.
% L. Baars       | 2023-Mar-01 | Added isCronJob feature, improved logging
%                |             | outputs, and added calculations of the HBR
%                |             | based on DISCOS reported shapes instead of
%                |             | just the span.
% L. Baars       | 2023-Mar-30 | Simplified the shape-based HBR
%                |             | calculations.
% L. Baars       | 2023-Mar-15 | Improved documentation for public release.
%                |             | Changed the isCronJob parameter to
%                |             | excludeBackspaceInLogOutputs.
% L. Baars       | 2025-Aug-04 | Added a filter to remove objects within
%                |             | reserved object ID ranges.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

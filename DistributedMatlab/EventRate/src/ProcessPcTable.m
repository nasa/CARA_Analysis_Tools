function params = ProcessPcTable(params)
% ProcessPcTable - Do the processing required to use the EventRate PcTable 
%                  operations mode
%
% Syntax: params = ProcessPcTable(params)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Do the processing required to use the EventRate PcTable 
%              operations mode
%
% =========================================================================
%
% Input:
%
%   params - EventRate parameter structure. See EventRate_default_params.m
%            or EventRate_ConjDist_default_params.m for documentation 
%
% =========================================================================
%
% Output:
%
%   params - EventRate parameter structure augmented with PcTable and OCMDB
%            information, depending on which mode EventRate is running in
%
% =========================================================================
%
% Dependencies:
%
%   BuildPcTable.m (if PcTable file is not provided)
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'../data')); addpath(s.path);
    s = what(fullfile(p,'../../../DistributedMatlab/Utils/LoggingAndStringReporting')); addpath(s.path);
    pathsAdded = true;
end

disp(' ');
disp([current_timestring() ' Processing OCMDB Pc table file']);

% Ensure a valid PcTable_mode has been specified
if min(abs(params.PcTable_mode-[0 1 2 3])) > 0
    error('Invalid PcTable_mode parameter');
end

% Check if we have been provided a PCTable file without an OCMDB file
if isempty(params.OCMDBfile) 
    % If we don't have an OCMDB file, look for a PcTable file
    if ~isempty(params.PcTableFile)
        NeedToBuildPcTable = false;
        params.haveOCMDBfile = false;
        PT = load(params.PcTableFile);
    else
        error('Error: Must at least provide either an OCMDB file or a PcTable file. \n')
    end
else
    params.haveOCMDBFile = true;
    % Check if OCMDB file name specifies the Pc table file or the original
    % OCMDB file
    PcTableString = '_PcTable';
    k = strfind(params.OCMDBfile,PcTableString);
    if isempty(k)
        % OCMDBfile specifies OCMDB file rather than Pc table file
        [p,f,e] = fileparts(params.OCMDBfile);
        f = [f PcTableString];
        params.PcTableFile = fullfile(p,[f e]);
    else
        % OCMDBfile specifies Pc table file rather than OCMDB file
        params.PcTableFile = params.OCMDBfile;
        params.OCMDBfile = strrep(params.PcTableFile,PcTableString,'');
    end
    % Check if Pc table file needs to be built
    if exist(params.PcTableFile,'file') 
        disp([current_timestring() ' Found OCMDB Pc table file: ' [f e]]);
        if params.PcTable_mode == 1
            % Load the old Pc table
            disp(' Loading OCMDB Pc table file....');
            PT = load(params.PcTableFile);
            if PT.PcTable.HBRmin_meters == params.PcTable_HBRmin_meters && ...
               PT.PcTable.HBRmax_meters == params.PcTable_HBRmax_meters && ...
               PT.PcTable.HBRnum        == params.PcTable_HBRnum
                NeedToBuildPcTable = false;
            else
                NeedToBuildPcTable = true;
                disp(' HBR table parameter mismatch - rebuilding Pc table file');
            end
        elseif params.PcTable_mode == 2
            NeedToBuildPcTable = true;
            disp(' Not using the previously built Pc table file');
        else
            error('Invalid PcTable_mode parameter');
        end
    else
        NeedToBuildPcTable = true;
    end
end


if NeedToBuildPcTable
    % Build the Pc table file, if required
    [params.PcTable,params.DB] = ...
        BuildPcTable(params.OCMDBfile,params.PcTableFile,params);
else
    % Copy the PcTable information from the PT structure
    % to the params structure
    params.PcTable = PT.PcTable;
    % Create the DB array in the params structure, replacing the original
    % deleted columns with NaN values
    Nrow = size(PT.DBTable,1);
    params.DB = NaN([Nrow PT.PcTable.ExpandedDBColumns]);
    kept = true([1 PT.PcTable.ExpandedDBColumns]);
    kept(PT.PcTable.DeletedDBColumns) = false;
    params.DB(:,kept) = PT.DBTable;
end

return;
end
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
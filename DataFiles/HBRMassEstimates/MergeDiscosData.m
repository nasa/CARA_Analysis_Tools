function MergeDiscosData(hbrMassFile)
% MergeDiscosData - Downloads and merges DISCOS satellite size and mass
%                   information into a NASA CARA HBR and Mass estimates
%                   file
%
% Syntax: MegeDiscosData(hbrMassFile);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   This method will connect to the European Space Agency (ESA) Database
%   and Information System Characterising Objects in Space (DISCOS) website
%   to download satellite size and mass information. It will then merge the
%   DISCOS data into the NASA CARA HBR and Mass information table for any
%   HBR or Mass entries labeled with 'DISCOS DB' data source labels. The
%   process will then create an output CSV file which includes the updated
%   HBR and Mass estimates.
%
%   The download and merge process is required to be performed by end users
%   since NASA CARA does not have the licensing data rights to redistribute
%   DISCOS data directly. This process is intended as a helper tool for end
%   users so that they don't have to create their own download and merge
%   processes.
%
% =========================================================================
%
% Input:
%
%    hbrMassFile - Path and file name to the CSV file containing the "for
%                  distribution" version of HBR and Mass estimates created
%                  by the NASA CARA team. The file name should have the
%                  format:
%                    HBR_Mass_Estimates_yyyymmdd_ForDistribution.csv
%
%                  The 'yyyymmdd' is a date tag associated with the HBR and
%                  mass estimate data.
%
%                  Upon successful completion of the script, an output file
%                  with merged DISCOS and NASA CARA data will be placed in
%                  the same directory as the original file, with a file
%                  name of:
%                    HBR_Mass_Estimates_yyyymmdd_Merged_NotForDistribution.csv
%
% =========================================================================
%
% Output:
%
%    None - see the description of the input parameter for the description
%           of the output file that will be created
%
% =========================================================================
%
% References:
%
%   Baars, L., and Hall, D., "Processing Space Fence Radar Cross-Section
%   Data to Produce Size and Mass Estimates," AAS Astrodynamics Specialist
%   Conference, AAS Paper 22-586, Aug. 2022.
%
%   Hall, D., and Baars, L., "Satellite Collision and Fragmentation
%   Probabilities Using Radar-Based Size and Mass Estimates," Journal of
%   Spacecraft and Rockets, Volume 60, Issue 4, Jul. 2023, pg. 1319-1332.
%
% =========================================================================
%
% Initial version: Nov 2023;  Latest update: Oct 2024
%
% ----------------- BEGIN CODE -----------------

    % Add pathing to the dummy password manager
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p,'DummyPasswordManager')); addpath(s.path);
        pathsAdded = true;
    end

    % Check the input filename to verify the proper format and create an
    % output filename
    [fPath, fName, fExt] = fileparts(hbrMassFile);
    if ~startsWith(fName,'HBR_Mass_Estimates_','IgnoreCase',true) || ...
            ~endsWith(fName,'_ForDistribution','IgnoreCase',true) || ...
            ~strcmpi(fExt,'.csv') || ...
            length(fName) ~= 43 || ...
            regexp(fName,'\d{8}') ~= 20
        error('This utility can only be used to merge DISCOS data into files named ''HBR_Mass_Estimates_yyyymmdd_ForDistribution.csv''');
    end
    newFName = strrep(fName,'_ForDistribution','_Merged_NotForDistribution');
    newHbrMassFileName = fullfile(fPath,[newFName fExt]);
    
    % Extract the date string from the filename and create the name of the
    % file which will store DISCOS data
    hbrMassDateStr = hbrMassFile(end-27:end-20);
    discosFile = fullfile(fPath,['discos_data_' hbrMassDateStr '.csv']);
    
    % Get the connection parameters
    connParams = connect_params_distrib;
    
    % Check for default Token
    secKey = GetSecKey;
    discosToken = GetPassword('DISCOS',secKey,connParams);
    if strcmp(discosToken,'<replace with token>')
        errMsg = ['Default token detected!' newline ...
            'Please follow the instructions in the "Setup Instructions" section' newline ...
            'of 01_Distribution_README.txt in order to setup a DISCOSweb access' newline ...
            'token for this tool.'];
        error(errMsg);
    end
    
    % Start the timing run and download the DISCOS data, store the data
    % into the DISCOS filename passed in
    tic
    discosData = GetDiscosData(discosFile, connParams);
    numDiscosRows = height(discosData);
    % Make sure the DISCOS data is sorted in ascending order by object ID
    discosData = sortrows(discosData,'ObjectID');
    
    % Read the input HBR/Mass file
    hbrMassData = readtable(hbrMassFile);
    % Make sure the HBR/Mass data is sorted in ascending order by object ID
    hbrMassData = sortrows(hbrMassData,'ObjectID');
    % Set some variables for hardcoded values
    discosID = 'DISCOS DB';
    defaultTxt = 'Default';
    % Setting to the defaults for unknown objects
    defaultHBR = 1.1;  % meters
    defaultMass = 448.3;  % kg
    
    % Find the indexes where either the HBR data source or the mass data
    % source indicate that the DISCOS DB should be the data source
    updateHBR = strcmpi(hbrMassData.HBRSource,discosID);
    updateMass = strcmpi(hbrMassData.MassSource,discosID);
    % Set a counter to keep track of the current index within the DISCOS
    % data
    j = 1;
    % Set counters for number of merged data points and number of data
    % points overridden with default values
    mergedHBR = 0;
    mergedMass = 0;
    setDefaultHBR = 0;
    setDefaultMass = 0;
    % Assuming that both data sources are sorted in ascending order based
    % on object ID, find the associated DISCOS entries and put their data
    % into the HBR/Mass data. This algorithm traverses both arrays in
    % parallel and is efficient since each array is traversed only once.
    % Loop through all of the HBR/Mass data
    for i = 1:height(hbrMassData)
        % If the HBR or mass needs to be updated with DISCOS data
        if updateHBR(i) || updateMass(i)
            % Get the current HBR/Mass object ID
            objectID = hbrMassData.ObjectID(i);
            % Traverse the DISCOS data from it's current position until the
            % DISCOS object ID is >= HBR/Mass object ID or until we've
            % gotten to the end of the DISCOS table
            while discosData.ObjectID(j) < objectID && j < numDiscosRows
                j = j + 1;
            end
            % If the DISCOS object ID matches the HBR/Mass object ID
            if discosData.ObjectID(j) == objectID
                % Merge in the DISCOS data for an HBR update
                if updateHBR(i)
                    if ~isnan(discosData.HBR(j))
                        hbrMassData.HBR(i) = discosData.HBR(j);
                        mergedHBR = mergedHBR + 1;
                    % In the unlikely event that the DISCOS HBR is NaN,
                    % override with the default HBR value
                    else
                        hbrMassData.HBRSource{i} = defaultTxt;
                        hbrMassData.HBR(i) = defaultHBR;
                        setDefaultHBR = setDefaultHBR + 1;
                    end
                end
                % merge in the DISCOS data for a Mass udpate
                if updateMass(i)
                    if ~isnan(discosData.Mass(j))
                        hbrMassData.Mass(i) = discosData.Mass(j);
                        mergedMass = mergedMass + 1;
                    % In the unlikely event that the DISCOS Mass is NaN,
                    % override with the default Mass value
                    else
                        hbrMassData.MassSource{i} = defaultTxt;
                        hbrMassData.Mass(i) = defaultMass;
                        setDefaultMass = setDefaultMass + 1;
                    end
                end
            % The object ID doesn't match, therefore the object ID doesn't
            % exist in the DISCOS table
            else
                % Override the HBR and mass values for DISCOS marked
                % updates if the object ID has not been found in the DISCOS
                % table
                if updateHBR(i)
                    hbrMassData.HBRSource{i} = defaultTxt;
                    hbrMassData.HBR(i) = defaultHBR;
                    setDefaultHBR = setDefaultHBR + 1;
                end
                if updateMass(i)
                    hbrMassData.MassSource{i} = defaultTxt;
                    hbrMassData.Mass(i) = defaultMass;
                    setDefaultMass = setDefaultMass + 1;
                end
            end
        end
    end
    
    % Display timing information and provide some information on the number
    % of values merged or overridden
    toc
    disp(['Merged ' num2str(mergedHBR) ' HBR values from the DISCOS DB']);
    disp(['Merged ' num2str(mergedMass) ' mass values from the DISCOS DB']);
    disp(['Set ' num2str(setDefaultHBR) ' default HBR values due to missing or NaN DISCOS DB values']);
    disp(['Set ' num2str(setDefaultMass) ' default mass values due to missing or NaN DISCOS DB values']);
    writetable(hbrMassData,newHbrMassFileName);
    disp(['Wrote merged data into the file ' newHbrMassFileName]);
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Nov-15 | Initial Development.
% L. Baars       | 2024-Oct-09 | Added a check for a default DISCOS token
%                                with better error handling. Also updated
%                                the default HBR/mass values to use the
%                                "unknown" values instead of "debris"
%                                values.

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

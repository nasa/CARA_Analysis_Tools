function [hbr, hbrSigma, hbrSource, mass, massSigma, massSource] = GetHBRMassEstimate(hbrMassFile, objectID)
% GetHBRMassEstimate - Gets size and mass information from the HBR and mass
%                      file
%
% Syntax: [hbr, hbrSigma, hbrSource, mass, massSigma, massSource] = GetHBRMassEstimate(hbrMassFile, objectID);
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
%   This method provides a convenient method for reading HBR and mass
%   information from an HBR and mass file. It has special checking in place
%   to indicate to a user that unmerged DISCOS data is used (causes the
%   function to exit with errors) or if the object does not have HBR/Mass
%   information or has invalid HBR/Mass information (both instances cause a
%   warning message and provide default HBR and mass values).
%
%   An end-user must pass in a merged HBR/Mass file that has been created
%   by the MergeDiscosData.m process.
%
% =========================================================================
%
% Inputs:
%
%   hbrMassFile - Name of the merged DISCOS HBR/Mass file.
%
%   objectID - NORAD catalog ID for the desired satellite
%
% =========================================================================
%
% Outputs:
%
%   hbr - Hard Body Radius (HBR) estimate of the satellite (in meters)
%   hbrSigma - Standard deviation of the HBR estimate (in meters). Values
%              are provided for hbrSource of 'SFK RCS Estimate', otherwise
%              this value will be 0.
%   hbrSource - Data source for the HBR estimate. Valid values are:
%                 CARA Protected Primary - HBR used for a NASA CARA mission
%                 DISCOS DB - HBR calculated from DISCOS size data
%                 SFK RCS Estimate - HBR estimated from Space Fence RCS
%                                    data
%                 Default - The value that is provided if no other data
%                           sources are found for the object
%   mass - Mass estimate of the satellite (in kg)
%   massSigma - Standard deviation of the mass estimate (in kg). Values are
%               provided for massSource of 'SFK RCS Estimate', otherwise
%               this value will be 0.
%   massSource - Data source for the Mass estimate. Valid values are the
%                same as the HBR source.
%   discosData - An nx14 table containing information about each object
%                retrieved from the DISCOS database. Field included are:
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

    % Read in the HBR and mass table, but use persistent variables so that
    % it only needs to occur once if the method is called again with the
    % same hbrMassFile
    persistent loadedFile hbrMassTable
    discosSrcTxt = 'DISCOS DB';
    if isempty(loadedFile) || ~strcmp(loadedFile, hbrMassFile)
        hbrMassTable = readtable(hbrMassFile);
        % Make sure the table has all the columns that we expect it to have
        hbrCols = hbrMassTable.Properties.VariableNames;
        if sum(ismember(hbrCols,'ObjectID')) == 0 || ...
                sum(ismember(hbrCols,'HBR')) == 0 || ...
                sum(ismember(hbrCols,'HBRSigma')) == 0 || ...
                sum(ismember(hbrCols,'HBRSource')) == 0 || ...
                sum(ismember(hbrCols,'Mass')) == 0 || ...
                sum(ismember(hbrCols,'MassSigma')) == 0 || ...
                sum(ismember(hbrCols,'MassSource')) == 0
            error(['Table defined in ' hbrMassFile ' is not a valid HBR/Mass table!']);
        end
        % Make sure that there are no DISCOS entries that haven't been
        % merged
        discosHBRIdxs = strcmp(hbrMassTable.HBRSource, discosSrcTxt);
        discosHBRs = hbrMassTable.HBR(discosHBRIdxs);
        discosMassIdxs = strcmp(hbrMassTable.MassSource, discosSrcTxt);
        discosMasses = hbrMassTable.Mass(discosMassIdxs);
        if sum(isnan(discosHBRs)) > 0 || sum(isnan(discosMasses)) > 0
            error(['Found unmerged DISCOS data within ' hbrMassFile '!' newline ...
                'Please run the MergeDiscosData.m utility to create a merged HBR/Mass estimates file' newline ...
                'and use that file as the input to this function.']);
        end
        loadedFile = hbrMassFile;
    end
    
    % Create the output Object ID string
    objectIDStr = sprintf('%09d',objectID);
    % Find the data for the object passed in
    tempTable = hbrMassTable(hbrMassTable.ObjectID == objectID,:);
    % Set variables for hardcoded values
    defaultTxt = 'Default';
    % Provide the defaults for "unknown" objects
    defaultHBR = 1.1;     % meters
    defaultMass = 448.3;  % kg
    % Provide a warning and return default values if the object isn't found
    if isempty(tempTable)
        warning(['No HBR/Mass estimates for object ' objectIDStr ' found in ' hbrMassFile newline ...
            'Providing default values for this object']);
        hbr = defaultHBR;
        hbrSigma = 0;
        hbrSource = defaultTxt;
        mass = defaultMass;
        massSigma = 0;
        massSource = defaultTxt;
        return;
    elseif height(tempTable) > 1
        % Provide a warning if multiple versions of the object has been
        % found in the table, only use the last instance of the object
        warning(['Multiple rows for object ' objectIDStr ' found in ' hbrMassFile newline ...
            'Using data from the last row in the table']);
        tempTable = tempTable(end,:);
    end
    
    % Provide warnings if the HBR or mass information are NaN, this should
    % never happen. Set these values to defaults.
    if isnan(tempTable.HBR)
        warning(['NaN HBR found for object ' objectIDStr ', setting the value to the default HBR (' num2str(defaultHBR) ' m)']);
        tempTable.HBR = defaultHBR;
        tempTable.HBRSigma = 0;
        tempTable.HBRSource{1} = defaultTxt;
    end
    if isnan(tempTable.Mass)
        warning(['NaN Mass found for object ' objectIDStr ', setting the value to the default Mass (' num2str(defaultMass) ' kg)']);
        tempTable.Mass = defaultMass;
        tempTable.MassSigma = 0;
        tempTable.MassSource{1} = defaultTxt;
    end
    % Provide warnings and set to 0 and sigma values that are NaN
    if isnan(tempTable.HBRSigma)
        warning(['NaN HBR sigma found for object ' objectIDStr ', setting the sigma to 0']);
        tempTable.HBRSigma = 0;
    end
    if isnan(tempTable.MassSigma)
        warning(['NaN Mass sigma found for object ' objectIDStr ', setting the sigma to 0']);
        tempTable.MassSigma = 0;
    end
    
    % Return the data that was found.
    hbr = tempTable.HBR;
    hbrSigma = tempTable.HBRSigma;
    hbrSource = tempTable.HBRSource{1};
    mass = tempTable.Mass;
    massSigma = tempTable.MassSigma;
    massSource = tempTable.MassSource{1};
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
% L. Baars       | 2024-Oct-09 | Updated the defaults to use the "unknown"
%                                type instead of the "debris" type.

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
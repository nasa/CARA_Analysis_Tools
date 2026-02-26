function CombinePartitionedFiles(inputDirName, outputDirName)
% CombinePartitionedFiles - Combines files named "<baseName>_Part*.mat"
%                           from the input directory into a single file
%                           named "<baseName>.mat" in the output directory,
%                           but only if the file doesn't already exist.
%
% Syntax: CombinePartitionedFiles(inputDirName);
%         CombinePartitionedFiles(inputDirName, outputDirName);
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
%   inputDirName - Name of the input directory which contains files named
%                  <baseName>_Part*.mat
%
%   outputDirName - (Optional) Name of the output directory where combined
%                   files will be written. All files with the same
%                   <baseName> will be combined into a file <baseName>.mat.
%                   (default = inputDirName)
%
% =========================================================================
%
% Output:
%
%   None - informational text is displayed to the console
%
% =========================================================================
%
% Initial version: Feb 2026;  Latest update: Feb 2026
%
% ----------------- BEGIN CODE -----------------

    % Check input parameters
    if nargin == 1
        outputDirName = inputDirName;
    elseif nargin ~= 2
        error('Incorrect number of parameters passed in');
    end

    % Check existence of input and output directories
    if ~exist(inputDirName,'dir')
        error(['Input directory not found: ' inputDirName]);
    end
    if ~exist(outputDirName,'dir')
        error(['Output directory not found: ' outputDirName]);
    end

    % Get all partitioned files matching the pattern "*_Part??.mat"
    fileNames = dir(fullfile(inputDirName,'*_Part*.mat'));
    numFiles = length(fileNames);
    if numFiles == 0
        warning(['Could not find any files named ' fullfile(inputDirName,'*_Part*.mat')]);
    end
    % Get the basenames for each of the files
    fileBaseNames = cell(numFiles,1);
    for i = 1:numFiles
        fileName = fileNames(i).name;
        fileBaseNames{i} = regexprep(fileName,'_Part\d\d\.mat$','');
    end
    
    % Load and combine the data for each unique basename
    uniqBaseNames = unique(fileBaseNames);
    numUniqBaseNames = length(uniqBaseNames);
    for i = 1:numUniqBaseNames
        outputFile = fullfile(outputDirName,[uniqBaseNames{i} '.mat']);
        if ~exist(outputFile,'file')
            disp(['Combining data for ' uniqBaseNames{i}]);
            allDBTable = [];
            fileNames = dir(fullfile(inputDirName,[uniqBaseNames{i} '_Part*.mat']));
            for j = 1:length(fileNames)
                fullFileName = fullfile(inputDirName, fileNames(j).name);
                if j == 1
                    % It is assumed that the first file will always include
                    % the PcTable
                    load(fullFileName,'DBTable','PcTable');
                else
                    load(fullFileName,'DBTable');
                end
                allDBTable = cat(1,allDBTable,DBTable);
                disp(['Loaded file ' fullFileName]);
            end
            DBTable = allDBTable;
            clear allDBTable;
            save(outputFile,'DBTable','PcTable','-v7.3');
            disp(['Created combined file: ' outputFile]);
        end
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% L. Baars  | 2026-Feb-24 | Initial version.
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
function BFMCMultiCA(parentDir)
% BFMCMultiCA - Multi-encounter Nc and Pc calculation using BFMC
%               relative distance minima tables.
%               (For CARA analysis team internal use)
%
% Syntax: BFMCMultiCA(parentDir);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Calculate Nc and Pc for a multi-encounter interaction using the
%   relative distance minima tables produced by the BFMC code (i.e.,
%   *.min files). Each rel. dist. minimum in the table represents a
%   separate conjunction, or close approach (CA) event, contained
%   within a separate encounter segment time period.
%
%   Specifically, the input files are from BFMC VCM-mode long-duration
%   runs.
%
% =========================================================================
%
% Input:
%
%    parentDir  -   Path to parent directory containing BFMC run output
%                   folders with *.min files
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

%% Initialization and defaults

% Add paths
persistent pathsAdded
if isempty(pathsAdded)
    addpath('..\..\DistributedMatlab\Utils\StateCovInterpolation');
    addpath('..\..\DistributedMatlab\Utils\AugmentedMath');
    addpath('..\..\DistributedMatlab\Utils\LoggingAndStringReporting');
    addpath('..\..\DistributedMatlab\Utils\Plotting');
    addpath('..\..\DistributedMatlab\Utils\CovarianceTransformations');
    addpath('..\..\DistributedMatlab\ProbabilityOfCollision');
    addpath('..\..\DistributedMatlab\ProbabilityOfCollision\Utils');
    addpath('..\..\DistributedMatlab\ProbabilityOfCollision\Pc3D_Hall_Utils');
    pathsAdded = true;
end    

% % Add paths
% persistent pathsAdded
% if isempty(pathsAdded)
%     addpath('src');
%     SDKpath = 'C:\Users\DHall\dth\Analysis\SDK\DistributedMatlab';
%     if exist(SDKpath,'dir') == 7
%         addpath(genpath(SDKpath));
%     else
%         SDKpath = 'D:\CARA Repository\Analysis\SDK\DistributedMatlab';
%         if exist(SDKpath,'dir') == 7
%             addpath(genpath(SDKpath));
%         else
%             error('SDKpath directory not found');
%         end
%     end
%     % SDKpath = 'C:\Users\DHall\dth\Analysis\SDK\ResearchCode\ProbabilityOfCollision';
%     % if exist(SDKpath,'dir') == 7
%     %     addpath(genpath(SDKpath));
%     % else
%     %     SDKpath = 'D:\CARA Repository\Analysis\SDK\ResearchCode\ProbabilityOfCollision';
%     %     if exist(SDKpath,'dir') == 7
%     %         addpath(genpath(SDKpath));
%     %     else
%     %         error('SDKpath directory not found');
%     %     end
%     % end
%     pathsAdded = true;
% end    

%% Set up to process the BFMC runs

% Create output directory in a parent directory assumed to contain
% a long-duration BFMC run populated a rel. dist. minima file (*.min)
outputDir = fullfile(parentDir,'MultiCA');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

% Find the conjunction ID directories output by the BFMC system
fileNames = dir(fullfile(parentDir,'*conj*','*.min'));
NfileNames = numel(fileNames);

%% Perform the 3D-Nc method processing

% Define and allocate output variables
out.outputDir = outputDir;
out.fileNames = fileNames;
out.Table = [];

% Process each BFMC directory
for i = 1:NfileNames
    
    % Extract conjunction ID
    [~, conjID, ~] = fileparts(fileNames(i).folder);
    disp(' ');
    disp(conjID);
    
    currData = BFMCMultiCAProcess(conjID,fileNames(i).folder,outputDir);
    out.Table = [out.Table; currData];
    
end

writetable(out.Table,fullfile(outputDir,'AnalysisResults.csv'));

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
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
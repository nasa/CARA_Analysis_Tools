function out = BFMCEphPcSequence(parentDir,EphSkip,UseMinEph,UseEqEls)
% BFMCEphPcSequence - Run EphemerisPc on BFMC output folders
%                     (For CARA analysis team internal use)
%
% Syntax: out = BFMCEphPcSequence(parentDir);
%         out = BFMCEphPcSequence(parentDir, EphSkip);
%         out = BFMCEphPcSequence(parentDir, EphSkip, UseMinEph);
%         out = BFMCEphPcSequence(parentDir, EphSkip, UseMinEph, UseEqEls);
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
%   Calculate Eph-Nc using the primary/secondary ephemeris ECI tables
%   produced by the BFMC code (i.e., *.eci files), optionally combined
%   with the relative distance minima tables produced by the BFMC code
%   (i.e., *.min).
%
%   Specifically, the input files are from BFMC VCM-mode long-duration
%   runs.
%
% =========================================================================
%
% Input:
%
%    parentDir  -   Path to parent directory containing BFMC run output 
%                   folders with *.eci and *.min files
%
%    EphSkip    -   Number of points to skip when processing *.eci
%                   ephemeris tables (optional, default = 0)
%
%    UseMinEph  -   Flag to augment *.eci ephemeris tables with
%                   rel. dist. minima points
%                   (optional, default = true)
%
%    UseEqEls   -   Flag to use equinoctial elements
%                   (optional, default = true)
%
% =========================================================================
%
% Output:
%
%   out         -   Output structure with fields:
%                     .outputDir   - Output directory path
%                     .fileNames   - Found conjunction file info
%                     .EphSkip     - EphSkip values used
%                     .UseMinEph   - UseMinEph flags used
%                     .UseEqEls    - UseEqEls flags used
%                     .BFMCEphPcOut - EphemerisPc results 
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

%% Initialization and defaults

% % Add paths
% persistent pathsAdded
% if isempty(pathsAdded)
%     addpath('..\Utils');    
%     addpath('..\AnalyzeScreenings\src');
%     addpath('..\..\DistributedMatlab\Utils\PosVelTransformations');
%     addpath('..\..\DistributedMatlab\Utils\CovarianceTransformations');
%     addpath('..\..\DistributedMatlab\Utils\StateCovInterpolation');
%     addpath('..\..\DistributedMatlab\Utils\AugmentedMath');
%     addpath('..\..\DistributedMatlab\Utils\LoggingAndStringReporting');
%     addpath('..\..\DistributedMatlab\Utils\General');
%     addpath('..\..\DistributedMatlab\Utils\Plotting');
%     addpath('..\..\DistributedMatlab\ProbabilityOfCollision\Utils');
%     addpath('..\..\DistributedMatlab\ProbabilityOfCollision\Pc3D_Hall_Utils');
%     pathsAdded = true;
% end    

Nargin = nargin;
if Nargin < 1
    error('Insufficient input');
end

% Number of points to skip when processing *.eci ephemeris tables
if Nargin < 2; EphSkip = []; end
if isempty(EphSkip); EphSkip = 0; end

% Augment *.eci ephemeris tables with rel. dist. minima points
if Nargin < 3; UseMinEph = []; end
if isempty(UseMinEph); UseMinEph = true(size(EphSkip)); end

% Ue equinoctial elements for curvilinear state estimations
if Nargin < 3; UseEqEls = []; end
if isempty(UseEqEls); UseEqEls = true(size(EphSkip)); end

%% Set up to process the BFMC runs

% Number of runs required
Nrun = numel(EphSkip);

% Ensure UseMinEph flags have correct dimensions
Nmin = numel(UseMinEph);
if Nrun ~= Nmin
    if Nrun > 1 && Nmin == 1
        UseMinEph = repmat(UseMinEph,size(EphSkip));
    else
        error('Dimensions of EphSkip and UseMinEph incompatible');
    end
end

% Ensure UseEphOffsetElements flags have correct dimensions
Nmin = numel(UseEqEls);
if Nrun ~= Nmin
    if Nrun > 1 && Nmin == 1
        UseEqEls = repmat(UseEqEls,size(EphSkip));
    else
        error('Dimensions of EphSkip and UseEqEls incompatible');
    end
end

% Create output directory in a parent directory assumed to contain
% a long-duration BFMC run populated with primary-secondary ephemeris files
% (*.eci) and a rel. dist. minima file (*.min)
outputDir = fullfile(parentDir,'BFMCEphPc');
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

% Find the conjunction ID directories output by the BFMC system
fileNames = dir(fullfile(parentDir,'*conj*','*.min'));
NfileNames = numel(fileNames);

%% Perform the sequential EphPc method processing

% Define and allocate output variables
out.outputDir = outputDir;
out.fileNames = fileNames;
out.EphSkip = EphSkip;
out.UseMinEph = UseMinEph;
out.UseEqEls = UseEqEls;
out.BFMCEphPcOut = cell(NfileNames,Nrun);

% Set up BFMC mode parameters for EphemerisPc function
EphemerisPcPars.BFMCmode = true;

% Process each BFMC directory
for i = 1:NfileNames
    
    % Extract conjunction ID
    [~, conjID, ~] = fileparts(fileNames(i).folder);
    disp(' ');
    disp(conjID);
    
    % Perform requested runs
    for j = 1:Nrun
        disp(' ');
        disp(['EphSkip = ' num2str(EphSkip(j)) ...
            ' UseMinEph = ' num2str(UseMinEph(j)) ...
            ' UseEqEls = ' num2str(UseEqEls(j))]);
        disp(' ');
        EphemerisPcPars.conjID = conjID;
        EphemerisPcPars.UseMinEph = UseMinEph(j);
        EphemerisPcPars.EphSkip = EphSkip(j);
        EphemerisPcPars.UseEqEls = UseEqEls(j);
        out.BFMCEphPcOut{i,j} = EphemerisPc( ...
            fileNames(i).folder,outputDir,EphemerisPcPars);
        % out.BFMCEphPcOut{i,j} = BFMCEphPc(conjID, ...
        %     fileNames(i).folder,outputDir, ...
        %     EphSkip(j),UseMinEph(j),UseEqEls(j));
    end
    
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
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
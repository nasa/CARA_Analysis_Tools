%% Add path
clear functions;
[p,~,~] = fileparts(mfilename('fullpath'));
s = what(fullfile(p,'../')); addpath(s.path);

%% Run High Relative Velocity LEO Example
fprintf('Running High Relative Velocity LEO Example\n');
inputPath = fullfile(p,'HighRelVelLEO/Inputs');
outputPath0 = fullfile(p,'HighRelVelLEO', ['Outputs_' datestr(now,'mm_dd_yyyy')]);
mkdir(outputPath0);
params.HBR = 20;
EphPcOut_HighRelVelLEO = EphemerisPc(inputPath,outputPath0,params);
fprintf('High Rel. Velocity LEO Run Outputs saved to: %s\n', outputPath0);

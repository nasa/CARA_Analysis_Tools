%% Run examples for the function EvaluateLightPollution

% - MMT data needs to be downloaded and placed in the subdirectories in 
%   'data/MMT' in order for the examples to run. 
%   Follow these instructions:
%     - Navigate to 'https://mmt9.ru/satellites/'.
%     - Navigate to the constellation subfolders in 'data/MMT/', which 
%       each contain a 'IDs_[Constellation].txt' file. The first 
%       line is an array of satellite IDs corresponding to all satellites 
%       used in the examples for that constellation. To retrieve these 
%       satellites, copy the entire line, paste it into the website’s 'ID' 
%       field, and then click 'Search'.
%     - For each satellite, download magnitude data files by selecting the 
%       'T' icon in the 'RCS' column to 'Download all tracks.' Files will 
%       be saved in the format 'satellite_[####].txt'.
%     - Place the files into the appropriate subfolder in 'data/MMT/'. 

%% Single-shell constellation examples

% Starlink 1st orbital shell constellation based on measured brightness
fprintf('\nRunning Starlink1stShell_MeasuredBrightness\n')
EvaluateLightPollution('params/Starlink1stShell_MeasuredBrightness');

% OneWeb Phase 1 constellation based on measured brightness
fprintf('\nRunning OneWebPhase1_MeasuredBrightness\n')
EvaluateLightPollution('params/OneWebPhase1_MeasuredBrightness');

% Iridium 2nd Gen. constellation based on measured brightness
fprintf('\nRunning Iridium2ndGen_MeasuredBrightness\n')
EvaluateLightPollution('params/Iridium2ndGen_MeasuredBrightness');

% OneWeb Phase 1 analysis using Iridium satellites as analog objects
% (to demonstrate analysis for proposed constellations with no orbiting
% satellites)
fprintf('\nRunning OneWebPhase1_IridiumAsAnalog\n')
EvaluateLightPollution('params/OneWebPhase1_IridiumAsAnalog');

%% Multi-shell constellation examples

% Starlink full 1st and 2nd gen. constellations based on measured
% brightness
fprintf('\nRunning Starlink1stGen_MeasuredBrightness\n')
EvaluateLightPollution('params/Starlink1stGen_MeasuredBrightness');
fprintf('\nRunning Starlink2ndGen_MeasuredBrightness\n')
EvaluateLightPollution('params/Starlink2ndGen_MeasuredBrightness');

% OneWeb Phase 2 constellation based on measured brightness
fprintf('\nRunning OneWebPhase2_MeasuredBrightness\n')
EvaluateLightPollution('params/OneWebPhase2_MeasuredBrightness');

%% Multi-shell constellation required brightness reduction example

% OneWeb Phase 2 constellation analysis showing that an average satellite
% brightness reduction of 1.7 magnitudes changes the evaluated light
% pollution level from "very high" to "medium"
fprintf('\nRunning OneWebPhase2_ReducedBrightness\n')
EvaluateLightPollution('params/OneWebPhase2_ReducedBrightness');

%% Multi-shell constellation comparison example

% Starlink constellations compared to OneWeb and Iridium constellations
fprintf('\nRunning Starlink1stGen_CompareConstellations\n')
EvaluateLightPollution('params/Starlink1stGen_CompareConstellations');
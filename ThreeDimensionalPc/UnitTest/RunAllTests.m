%% RunAllTests
% This code is intended to run all written Unit Tests for the Software
% Development Kit (SDK)
%
% T. Lechtenberg            Feb 2018

%% Set up File Paths
p = mfilename('fullpath');
[filepath,~,~] = fileparts(p);
cd([filepath filesep '..']);
FullPath=pwd;
addpath(genpath(FullPath));
cd(filepath);

% Turn Warnings off
warning off
%% Get all test suites
list=dir([ '**' filesep '*.m']);

Tests=[];
for i=1:length(list)
    if ~strcmpi(list(i).name,'RunAllTests.m') %%&& ~strcmpi(list(i).name,'Pc_MC_Kep2body_parallel_UnitTest.m')
        Tests = [Tests, matlab.unittest.TestSuite.fromFile([list(i).folder filesep list(i).name])];
    end
end

%% Run Full Test Suite
result = run(Tests);

%% Print Results to Command Line
T = table(result)
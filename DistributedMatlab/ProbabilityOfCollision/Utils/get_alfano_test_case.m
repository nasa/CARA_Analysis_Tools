function [tcdata] = get_alfano_test_case(casenum,datapath)

% Read Excel files with Alfano (2009) test-case data
%
% INPUT:

%   casenum     = Case number (1 through 12)
%   datapath    = Directory for XLS file data (default = pwd)
%
% OUTPUT:
%
%   tcdata      = Structure with test case data
%
% REFERENCE:
%
% S.Alfano, "Satellite Conjunction Monte Carlo Analysis" AAS 09-233 (2009).
%

% Set the default path

if nargin < 2; datapath = []; end
if isempty(datapath); datapath = pwd; end

% Define the XLS data file names

casestr = ['Case' num2str(casenum) ' data at '];
file_epoch = fullfile(datapath,[casestr 'Epoch Time.xls']);
file_tca   = fullfile(datapath,[casestr 'TCA.xls']);

% Read the XLS data for EPOCH

[~,~,raw] = xlsread(file_epoch);
tcdata.casenum = casenum;

% Define the auxilliary data (not tabulated in the Excel files)

switch casenum
    case 1
        tcdata.HBR = 15;
        tcdata.PC_MC_10e8 = 0.217467140;
        tcdata.desc = 'GEO, nonlinear relative motion';
    case 2
        tcdata.HBR = 4;
        tcdata.PC_MC_10e8 = 0.015736620;
        tcdata.desc = 'GEO, nonlinear relative motion';
    case 3
        tcdata.HBR = 15;
        tcdata.PC_MC_10e8 = 0.100846420;
        tcdata.desc = 'GEO, linear relative motion';
    case 4
        tcdata.HBR = 15;
        tcdata.PC_MC_10e8 = 0.073089530;
        tcdata.desc = 'GEO, nonlinear relative motion';
    case 5
        tcdata.HBR = 10;
        tcdata.PC_MC_10e8 = 0.044498913;
        tcdata.desc = 'LEO, linear relative motion';
    case 6
        tcdata.HBR = 10;
        tcdata.PC_MC_10e8 = 0.004300500;
        tcdata.desc = 'LEO, near-linear relative motion';
    case 7
        tcdata.HBR = 10;
        tcdata.PC_MC_10e8 = 0.000161462;
        tcdata.desc = 'LEO, nonlinear relative motion';
    case 8
        tcdata.HBR = 4;
        tcdata.PC_MC_10e8 = 0.035256080;
        tcdata.desc = 'MEO, nonlinear relative motion';
    case 9
        tcdata.HBR = 6;
        tcdata.PC_MC_10e8 = 0.365116060;
        tcdata.desc = 'HEO, nonlinear relative motion';
    case 10
        tcdata.HBR = 6;
        tcdata.PC_MC_10e8 = 0.362952470;
        tcdata.desc = 'HEO, nonlinear relative motion';
    case 11
        tcdata.HBR = 4;
        tcdata.PC_MC_10e8 = 0.003328530;
        tcdata.desc = 'LEO, leader-follower';
    case 12
        tcdata.HBR = 4;
        tcdata.PC_MC_10e8 = 0.002555950;
        tcdata.desc = 'LEO, identical orbits';
    otherwise
        error('Case number out of range (1 through 12)');
end

% Extract the Excel data

n = 5;
tcdata.R1e  = [raw{n,1}; raw{n,2}; raw{n,3}];

n = 9;
tcdata.V1e  = [raw{n,1}; raw{n,2}; raw{n,3}];

n=13;
tcdata.P1e  = [raw{n  ,1} raw{n  ,2} raw{n  ,3} raw{n  ,4} raw{n  ,5} raw{n  ,6};
               raw{n+1,1} raw{n+1,2} raw{n+1,3} raw{n+1,4} raw{n+1,5} raw{n+1,6};
               raw{n+2,1} raw{n+2,2} raw{n+2,3} raw{n+2,4} raw{n+2,5} raw{n+2,6};
               raw{n+3,1} raw{n+3,2} raw{n+3,3} raw{n+3,4} raw{n+3,5} raw{n+3,6};
               raw{n+4,1} raw{n+4,2} raw{n+4,3} raw{n+4,4} raw{n+4,5} raw{n+4,6};
               raw{n+5,1} raw{n+5,2} raw{n+5,3} raw{n+5,4} raw{n+5,5} raw{n+5,6}];

n = 23;
tcdata.R2e  = [raw{n,1}; raw{n,2}; raw{n,3}];

n = 27;
tcdata.V2e  = [raw{n,1}; raw{n,2}; raw{n,3}];

n=31;
tcdata.P2e  = [raw{n  ,1} raw{n  ,2} raw{n  ,3} raw{n  ,4} raw{n  ,5} raw{n  ,6};
               raw{n+1,1} raw{n+1,2} raw{n+1,3} raw{n+1,4} raw{n+1,5} raw{n+1,6};
               raw{n+2,1} raw{n+2,2} raw{n+2,3} raw{n+2,4} raw{n+2,5} raw{n+2,6};
               raw{n+3,1} raw{n+3,2} raw{n+3,3} raw{n+3,4} raw{n+3,5} raw{n+3,6};
               raw{n+4,1} raw{n+4,2} raw{n+4,3} raw{n+4,4} raw{n+4,5} raw{n+4,6};
               raw{n+5,1} raw{n+5,2} raw{n+5,3} raw{n+5,4} raw{n+5,5} raw{n+5,6}];
           
% Read the XLS data for TCA

[~,~,raw] = xlsread(file_tca);

n = 5;
tcdata.R1o  = [raw{n,1}; raw{n,2}; raw{n,3}];

n = 9;
tcdata.V1o  = [raw{n,1}; raw{n,2}; raw{n,3}];

n=13;
tcdata.P1o  = [raw{n  ,1} raw{n  ,2} raw{n  ,3} raw{n  ,4} raw{n  ,5} raw{n  ,6};
               raw{n+1,1} raw{n+1,2} raw{n+1,3} raw{n+1,4} raw{n+1,5} raw{n+1,6};
               raw{n+2,1} raw{n+2,2} raw{n+2,3} raw{n+2,4} raw{n+2,5} raw{n+2,6};
               raw{n+3,1} raw{n+3,2} raw{n+3,3} raw{n+3,4} raw{n+3,5} raw{n+3,6};
               raw{n+4,1} raw{n+4,2} raw{n+4,3} raw{n+4,4} raw{n+4,5} raw{n+4,6};
               raw{n+5,1} raw{n+5,2} raw{n+5,3} raw{n+5,4} raw{n+5,5} raw{n+5,6}];

n = 23;
tcdata.R2o  = [raw{n,1}; raw{n,2}; raw{n,3}];

n = 27;
tcdata.V2o  = [raw{n,1}; raw{n,2}; raw{n,3}];

n=31;
tcdata.P2o  = [raw{n  ,1} raw{n  ,2} raw{n  ,3} raw{n  ,4} raw{n  ,5} raw{n  ,6};
               raw{n+1,1} raw{n+1,2} raw{n+1,3} raw{n+1,4} raw{n+1,5} raw{n+1,6};
               raw{n+2,1} raw{n+2,2} raw{n+2,3} raw{n+2,4} raw{n+2,5} raw{n+2,6};
               raw{n+3,1} raw{n+3,2} raw{n+3,3} raw{n+3,4} raw{n+3,5} raw{n+3,6};
               raw{n+4,1} raw{n+4,2} raw{n+4,3} raw{n+4,4} raw{n+4,5} raw{n+4,6};
               raw{n+5,1} raw{n+5,2} raw{n+5,3} raw{n+5,4} raw{n+5,5} raw{n+5,6}];
           
% Read the pre-tabulated TCA in seconds after EPOCH

str = raw{1,1};
[p,~] = string_parts(str);

if ~strcmpi('epoch)',p{end})
    tcdata.tca0_from_epoch = NaN;
else
    tcdata.tca0_from_epoch = str2double(p{end-3});
    
end

if isnan(tcdata.tca0_from_epoch)
    warning('Unable to read tabulated TCA');
end

% Calculate the linear estimates for the TCA

r12 = tcdata.R1e-tcdata.R2e;
v12 = tcdata.V1e-tcdata.V2e;
tcdata.tca_from_epoch_lin = -(r12'*v12)/(v12'*v12); 

r12 = tcdata.R1o-tcdata.R2o;
v12 = tcdata.V1o-tcdata.V2o;
tcdata.tca_from_tca0_lin  = -(r12'*v12)/(v12'*v12); 

return
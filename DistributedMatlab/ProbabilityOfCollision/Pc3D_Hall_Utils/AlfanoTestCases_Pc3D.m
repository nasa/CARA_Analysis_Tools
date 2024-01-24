function AlfanoTestCases_Pc3D(CaseNumbers)

% Process one or more conjunction test cases from Alfano (2009) using the
% function Pc3D_Hall, which implements the Hall (2021) "3D-Nc" algorithm
% that calculates Pc estimates for curvilinear and extended-duration
% satellite interactions.

% REFERENCES:
%
% S.Alfano, "Satellite Conjunction Monte Carlo Analysis" AAS 09-233, 2009.
% (also referred to as "A09")
%
% D.Hall, "Expected Collision Rates for Tracked Satellites" Journal of
% Spacecraft and Rockets, Vol. 58, No. 3, May–June 2021, pp. 715–728.
% https://doi.org/10.2514/1.A34919 (also referred to as "H21")
%
% D.Hall, L.Baars, S.Casali, "A Multistep Probability of Collision
% Computational Algorithm" AAS 23-398, 2023. (also referred to as HBC23)
%

%% Initializations and defaults

% Add required paths
addpath(genpath('../..'));

% Defaults
if nargin < 1; CaseNumbers = []; end

% Check case numbers
if isempty(CaseNumbers); CaseNumbers = (1:12); end
BadCaseNumbers = CaseNumbers < 1 | CaseNumbers > 12 | ...
                 round(CaseNumbers) ~= CaseNumbers;
if any(BadCaseNumbers)
    warning('Excluding invalid case numbers');
    CaseNumbers = CaseNumbers(~BadCaseNumbers);
end
if numel(CaseNumbers) == 0
    warning('No valid case numbers found; no processing attempted');
    return;
end

% Pc display format
PcFmt = '%0.6e';

%% Get the Alfano (2009) data for the test case conjunctions

% Read tables of A09 at-TCA conjunction data
pars.case_list = CaseNumbers;
pars.data_path = 'Input_AlfanoTestCases_Pc3D';
pars.verbose = true;
conj = GetAlfanoTestCases(pars);
Nconj = numel(conj.case);

% Get the supplementary data copied from the A09 paper, and calculated
% using the CARA Simplified Dynamics Monte Carlo method
SupDat = readtable(fullfile(pars.data_path, ...
                   'Alfano09_CARA_Supplementary_Data.xlsx'));
SupDat = SupDat(CaseNumbers,:);

% Alfano (2002) quantities, extracted from the tables in the A09 paper
A09_FinalTime   = SupDat.A09_FinalTime_s; % A09 time interval bounds
A09_InstProb    = SupDat.A09_InstProb;    % A09 instantaneous probability
A09_PcMC1e8     = SupDat.A09_PcMC1e8;     % A09 Two-Body Monte Carlo (TBMC) Pc estimate
A09_PcLinear100 = SupDat.A09_PcLinear100; % A09 2D-Pc estimate, using 100 integration points

% CARA quantities, pre-calculated using the methods
% described in the H21 and HBC23 papers
CARA_PcMC       = SupDat.CARA_PcMC;       % CARA TBMC (or SDMC) Pc estimate
CARA_PcMCLo95   = SupDat.CARA_PcMCLo95;   % CARA lower bound of MC 95% conf. range
CARA_PcMCHi95   = SupDat.CARA_PcMCHi95;   % CARA upper bound of MC 95% conf. range
CARA_MCHits     = SupDat.CARA_MCHits;     % CARA number of MC hits
CARA_MCTrials   = SupDat.CARA_MCTrials;   % CARA number of MC trials

%% Process the conjunctions

% Allocate arrays
sz   = [numel(CaseNumbers) 1];
Case = reshape(CaseNumbers,sz);
HBR  = NaN(sz);
CARA_Pc2D = NaN(sz);
CARA_Pc3D = NaN(sz);

% Loop over conjunctions
for nc=1:Nconj
    
    % Hard-body radius
    HBR(nc) = conj.HBR(nc);
    
    % States and covariances
    r1 = conj.X1(1:3,nc)'; v1 = conj.X1(4:6,nc)'; C1 = conj.C1(:,:,nc);
    r2 = conj.X2(1:3,nc)'; v2 = conj.X2(4:6,nc)'; C2 = conj.C2(:,:,nc);
    
    % Bounds for 3D-Nc method Pc calculation, to account for the A09
    % temporal encounter segments that are very long in duration
    TLimit                  = A09_FinalTime(nc);
    Pc3DParams.Tmin_initial = -TLimit;
    Pc3DParams.Tmax_initial =  TLimit;
    Pc3DParams.Tmin_limit   = Pc3DParams.Tmin_initial;
    Pc3DParams.Tmax_limit   = Pc3DParams.Tmax_initial;
    
    % Initial number of time integration points for 3D-Nc calculation,
    % also to account for the A09 long-duration encounter segments
    Pc3DParams.Neph = 1000;
    
    % Calculate the 2D-Pc estimate of the Pc value, which can be compared
    % to the A09_PcLinear100 value
    vmag = norm(v2-v1);
    if vmag ~= 0
        CARA_Pc2D(nc) = Pc2D_Foster(r1,v1,C1,r2,v2,C2,HBR(nc), ...
                                    1e-9,'circle');
        Pc2DMethod = '2D-Pc method using CARA Pc2D_Foster function';
    else
        Pc2DMethod = ...
            '2D-Pc method undefined for zero relative velocity conjunctions';
    end
    
    % Calculate the 3D-Nc method estimate of the Pc value
    [CARA_Pc3D(nc), Pc3DInfo] = Pc3D_Hall(r1,v1,C1,r2,v2,C2, ...
                                         HBR(nc),Pc3DParams);

    Pc3DMethod = ['3D-Nc method using CARA Pc3D_Hall function' ...
        ' with ' num2str(Pc3DInfo.Neph) ' time integration points'];

    % Display results
    disp(' ');
    disp(['==== Alfano (2009) case number ' num2str(conj.case(nc)) ' ====']);
    disp([' HBR (m) = ' num2str(HBR(nc))]);
    disp([' 2D-Pc = ' num2str(CARA_Pc2D(nc),PcFmt)]);
    disp(['   (' Pc2DMethod ')']);
    disp([' Time bounds (+/- s) = ' num2str(TLimit)]);
    disp([' 3D-Pc = ' num2str(CARA_Pc3D(nc),PcFmt)]);
    disp(['   (' Pc3DMethod  ')']);
    
    % Plot the cumulative Pc over the A09 encounter time segment
    figure(nc); clf;
    plot(Pc3DInfo.Teph,Pc3DInfo.Nccum,'.-');
    xlabel('t - TCA (s)');
    ylabel('Cumulative Pc');
    titl = {['Alfano (2009) case number ' num2str(conj.case(nc))], ...
            ['CARA 2D-Pc = ' num2str(CARA_Pc2D(nc),PcFmt)], ...
            ['CARA 3D-Pc = ' num2str(CARA_Pc3D(nc),PcFmt)], ...
            ' '};
    title(titl);
    drawnow;
    
end

%% Output tabulated results

% Check for output directory
OutDir = 'Output_AlfanoTestCases_Pc3D';
chkoutput=exist(OutDir,'dir');
if chkoutput == 0
    mkdir(OutDir);
elseif chkoutput ~= 7
    error(['Problem with output directory existence: ' OutDir]);
end

% Create and write output tables in both .xlsx and .csv format
t = table(Case,HBR, ...
          A09_InstProb,A09_PcMC1e8,A09_PcLinear100, ...
          CARA_Pc2D,CARA_Pc3D,CARA_PcMC, ...
          CARA_PcMCLo95,CARA_PcMCHi95, ...
          CARA_MCHits,CARA_MCTrials);
f = ['A09TestCases_' num2str(min(CaseNumbers),'%02i') ...
    '_' num2str(max(CaseNumbers),'%02i')];
writetable(t,fullfile(OutDir,[f '.xlsx']));
writetable(t,fullfile(OutDir,[f '.csv']));

return
end
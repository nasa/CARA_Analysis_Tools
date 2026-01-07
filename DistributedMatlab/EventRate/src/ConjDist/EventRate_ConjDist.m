function [out,params] = EventRate_ConjDist(params)
% EventRate_ConjDist - Calculate and plot conjunction distributions for a 
%                      specified set of CARA primary objects
% Syntax: [out, params] = EventRate_ConjDist(params);
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculate and plot conjunction distributions for a specified 
% set of CARA primary objects
%
% This code has been adapted from the original ConjDist.m R&D code,
% adjusted to work with the EventRate.m function, by turning several
% of its original modes off.
%
% =========================================================================
%
% Input:
%
%   params - Auxilliary input parameter structure
%                       
% =========================================================================
%
% Output:
%
%   out    - Output structure with the following fields:
%
%               MissionEventRate: 5-50-95 percentile estimated red  
%               event rates for the given mission at the input Pc 
%               threshold 
%
%               ConjDist: Output substructure containing various 
%               intermediate values used to obtain output plots
%
%   params - Structure. See 'EventRate_ConjDist_default_params.m' for more
%            documentation
%
% =========================================================================
%
% References:
%
%   Doyle Hall (2019) "Determining Appropriate Risk Remediation Thresholds
%   from Empirical Conjunction Data using Survival Probability Methods" 
%   AAS 19-631.
%
% =========================================================================
%
% Dependencies:
%
%   ConjDist_binofit.m
%   ConjDist_set_resample.m
%   EventRate_ConjDist_default_params.m
%   ConjDist_update_plot.m
%   ConjDist_altlat_plot.m
%   ConjDist_time_plot.m (Only if reporting SVI rates)
%
%   ../../../CovarianceRealism
%   ../../../Utils/OrbitTransformations
%   ../../../Utils/CovarianceTransformations
%   ../../../Utils/LoggingAndStringReporting
%   ../../../Utils/Plotting
%   ../../../Utils/TimeTransformations
%
%   For 'exclude non-catastrophic events':
%
%   ../../../ProbabilityOfCollision
%   ../../../Utils/AugmentedMath
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Oct 2025
%
% ----------------- BEGIN CODE -----------------
%% Initializations and defaults

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'../matlab_lib')); addpath(s.path);
    s = what(fullfile(p,'../../../CovarianceRealism')); addpath(s.path);
    s = what(fullfile(p,'../../../CollisionConsequence')); addpath(s.path);
    s = what(fullfile(p,'../../../Utils/OrbitTransformations')); addpath(s.path);
    s = what(fullfile(p,'../../../Utils/CovarianceTransformations')); addpath(s.path);
    s = what(fullfile(p,'../../../Utils/LoggingAndStringReporting')); addpath(s.path);
    s = what(fullfile(p,'../../../Utils/Plotting')); addpath(s.path);
    s = what(fullfile(p,'../../../Utils/TimeTransformations')); addpath(s.path);
    pathsAdded = true;
end

Nargin = nargin;
if (Nargin < 1); params =[]; end

% Flag indicating if pretabulated Pc table mode is being used
using_PcTable = ~isempty(params.PcTable);

% Matlab date number zero point JD
JD0 = timestring_to_jd('0000-01-00 00:00:00'); 

% Figure handles
fighandle0 = [];
fighandle1 = [];
fighandle2 = [];

% Set up for plotting
xfsz = 12;
xfwt = 'bold';
yfsz = xfsz;
yfwt = xfwt;
tfsz = xfsz;
tfwt = 'bold';
tfan = 'italic';
afsz = xfsz;
afwt = 'bold';
alwd = 1;

% Padding string
pad20 = repmat(' ',[1 20]);

% 95% confidence chi-squared factor
chifact = sqrt(chi2inv(0.95,1));

% Initialize output
out = []; 

% Initialize PcMultiStep parameters
PcMultiStep_params = [];
PcMultiStep_params.OnlyPc2DCalculation = true;

% Set up EventRate processing modes

if params.EventRate_estimation_mode == 1
    % Given Tmis and Prem, find Pcum(noRMM), Pcum(w/RMM) and RMM-Rate
    PcRemEstimation = false;
    params.Pcmission = NaN;
    params.PcRMMexecution = params.EventRate_Pc_value;
    L10PremExecution = log10(params.PcRMMexecution);
elseif params.EventRate_estimation_mode == 2
    % Given Tmis & Pcum(w/RMM), find Pcum(noRMM) and Prem
    PcRemEstimation = true;
    params.Pcmission = params.EventRate_Pc_value;
    params.PcRMMexecution = NaN;
    L10Pcmission = log10(params.Pcmission);
else
    error('Invalid EventRate_estimation_mode parameter');
end

% Other ConjDist parameters
params.Tmissionyrs = params.mission_duration_years;
params.CommitConsiderHours = 24*[params.Tcommit_days params.Tconsider_days];
params.add_one_unrem_event = false;

% Prem image array
L10Pr1 = -7; L10Pr2 = -2;
if ~PcRemEstimation
    LPcRMMexecution = log10(params.PcRMMexecution);
    L10Pr1 = min(L10Pr1,LPcRMMexecution);
    L10Pr2 = max(L10Pr2,LPcRMMexecution);
end
NL10Pr = round(10*(L10Pr2-L10Pr1))+1;
params.L10PremImage = [L10Pr1 L10Pr2 NL10Pr];
params.log10PremSearch = params.L10PremImage;

% Make Rho vs Prem image (only in mode 2 when estimating Prem)
params.FredVsPremImage = true;
params.L10FredImage = [-2.5 0 26]; % [L10Fr1 L10Fr2 NFr]

% Table of mission durations
params.TmisTable = [0.5 1 2 3 5 7.5 10 12.5 15 20 25];
params.L10PcmisTable = [-6 -5 -4 -3 -2];

% Parameters of mission duration array for images (must be constructed to
% include all TmisTable entries)
params.TmisImage = [0.5 25 50]; % [Tm1 Tm2 NTm]

% Parameters for mission duration array for x-y plots
params.Tmissionyrsplot = [params.TmisImage(1)/1.5 params.TmisImage(2)*1.5 params.TmisImage(3)*2];

% Make conjunction latitude vs altitude spatial distribution plots
if params.plot_spatial_distribution
    params.make_altlat_plots = -4;
else
    params.make_altlat_plots = 0;
end

if strcmpi(params.redyel_event_option,'lud')
    % Last update Pc value between consider and commit times, 
    % used for modeling RMM execution rates
    PcLUDvsPcMaxOption = 1;
elseif strcmpi(params.redyel_event_option,'mud')
    % Maximum update Pc value between consider and commit times, 
    % used for modeling RMM planning rates
    PcLUDvsPcMaxOption = 2;
else
    % Non standard operational mode
    PcLUDvsPcMaxOption = 0;
end

% Construct an output file root name
if PcRemEstimation
    ModeStr = 'RMMPcLevel';
else
    if PcLUDvsPcMaxOption == 1
        ModeStr = 'RMMExecutionRate';
    elseif PcLUDvsPcMaxOption == 2
        ModeStr = 'RMMPlanningRate';
    else
        ModeStr = 'EventRate';
    end
end

MisNamStr = ['_' params.mission_name];

if isfield(params,'priorb')
    if ~isempty(params.priorb)
        OrbStr = ['_' params.priorb];
    end
else
    OrbStr = '';
end

% Base output file root
output_file_root0 = [ModeStr MisNamStr OrbStr '_' params.prisetstr ...
    '_HBR' num2str(params.mission_HBR_meters) 'm' ...
    '_Com' smart_exp_format(params.Tcommit_days,3) 'to' ...
           smart_exp_format(params.Tconsider_days,3) 'dy'];

% If excluding likely non-catastrophic collisions, include mass
% estimate in output file names
if params.exclude_noncatastrophic
    output_file_root0 = [output_file_root0 '_Mass' ...
        smart_exp_format(params.mission_on_orbit_mass_kg,3) 'kg'];
end

% Add mission duration to output file names
output_file_root = [output_file_root0 ...
                  '_Dur' num2str(params.Tmissionyrs) 'yr'];

% Add secondary catalog growth parameter to output file names
SecCatGrowth = params.SecondaryCatalogGrowth;
output_file_root = [output_file_root '_Grow' num2str(SecCatGrowth)];

% Add conservative mode flag to output file names
if params.add_one_unrem_event
    output_file_root = [output_file_root '_Conservative'];
end

% Set the run tag to be the same as the output file root
params.runtag = output_file_root;

% Set remaining ConjDist parameters to defaults
params = EventRate_ConjDist_default_params(params);

% Turn off EventRate_ConjDist flags that are incompatible with PcTable mode
if using_PcTable
    if params.association_check
        warning('PcTable incompatible with association_check parameter - turning off');
        params.association_check = false;
    end
end

% Log/display initiation time
logstr = ' ';
log_string(params.logfid,logstr,params.logging,params.displaying);
logstr = [params.start_time ' Executing function EventRate_ConjDist'];
log_string(params.logfid,logstr,params.logging,params.displaying);

%% Load database (DB) contained in the OCMDB file, if required

if using_PcTable
    % If we're in PcTable mode, we don't have to check for an OCMDB file
    logstr = [current_timestring() ' Using OCMDB Pc table processing mode'];
    log_string(params.logfid,logstr,params.logging,params.displaying);
else
    if isempty(params.DB)
        logstr = [current_timestring() ' Loading OCMDB file ' ...
            params.OCMDBroot params.OCMDBext];
        log_string(params.logfid,logstr,params.logging,params.displaying);
        params.DB = load(params.OCMDBfile);
        params.DB = params.DB.DB;
    else
        logstr = [current_timestring() ' Found DB array for OCMDB file ' ...
            params.OCMDBroot params.OCMDBext];
        log_string(params.logfid,logstr,params.logging,params.displaying);
    end
end

[Ntot,Ncol] = size(params.DB);
if (Ncol < 218)
    error('Insufficient number of OCMDB columns');
end
logstr = [current_timestring() ' Conjunctions in the entire DB: ' ...
    num2str(Ntot)];
log_string(params.logfid,logstr,params.logging,params.displaying);


if ~using_PcTable
    % If not using a PcTable, then eliminate bad conjunctions.
    % (This elimination has already been done if using a PcTable.)
    params.DB = EliminateBadDBEntries(params.DB,1,1,1,0);
    Ntot = size(params.DB,1);
    logstr = [current_timestring() ...
        ' Conjunctions in DB after eliminating bad entries: ' ...
        num2str(Ntot)];
    log_string(params.logfid,logstr,params.logging,params.displaying);
end

% Report the primary set
logstr = [current_timestring() ' Primary set ID = ' params.prisetstr];
log_string(params.logfid,logstr,params.logging,params.displaying);
Npri = params.Npri;
if Npri > 1
    for np=1:Npri
        logstr = [pad20 num2str(params.priset(np),'%05i ')];
        log_string(params.logfid,logstr,params.logging,params.displaying);
    end
end

% Copy the DB over to an output variable before adjustments
out.ConjDist.DB = params.DB;

% Build the covariance matrix blocks
CovBuild6 = triu(ones(6,6))';
CovBuild6(CovBuild6~=0)  = 1:21;
CovBuild6 = CovBuild6 + triu(CovBuild6',1);

% Record primary set in output
out.ConjDist.Npri = Npri;
out.ConjDist.priset = params.priset; out.ConjDist.prisetstr = params.prisetstr;

% Use DB TCA extrema as the TCA limits if not provided
if isempty(params.TCAlimits_UT)
    impose_TCAlimits = false;
    out.ConjDist.TCAlimits(1) = min(out.ConjDist.DB(:,166));
    out.ConjDist.TCAlimits(2) = max(out.ConjDist.DB(:,166));
else
    impose_TCAlimits = true;
    if strcmpi(class(params.TCAlimits_UT),'double')
        % Relative limits
        out.ConjDist.TCAlimits(1) = min(out.ConjDist.DB(:,166))+params.TCAlimits_UT(1);
        out.ConjDist.TCAlimits(2) = max(out.ConjDist.DB(:,166))-params.TCAlimits_UT(2);
    else
        % Absolute limits
        out.ConjDist.TCAlimits(1) = datenum(round_timestring(params.TCAlimits_UT{1}), ...
            'yyyy-mm-dd HH:MM:SS.FFF');
        out.ConjDist.TCAlimits(1) = max(out.ConjDist.TCAlimits(1),min(out.ConjDist.DB(:,166)));
        out.ConjDist.TCAlimits(2) = datenum(round_timestring(params.TCAlimits_UT{2}), ...
            'yyyy-mm-dd HH:MM:SS.FFF');
        out.ConjDist.TCAlimits(2) = min(out.ConjDist.TCAlimits(2),max(out.ConjDist.DB(:,166)));
        out.ConjDist.TCAlimits = sort(out.ConjDist.TCAlimits);
    end
end

% Eliminate all 90000 series secondaries
ndx = (out.ConjDist.DB(:,2) >= 90000);
if any(ndx)
    logstr = [current_timestring() ...
        ' Eliminating ' num2str(sum(ndx)) ...
        ' conjunctions involving 90000-series secondaries'];
    log_string(params.logfid,logstr,params.logging,params.displaying);    
    out.ConjDist.DB = out.ConjDist.DB(~ndx,:); 
end

% Find all conjunctions in the entire DB involving the specified primaries
% and secondaries

ndx = ismember(out.ConjDist.DB(:,1),params.priset);

if ~isempty(params.secset)
    secstr = '+secondary';
    ndx = ndx & ismember(out.ConjDist.DB(:,2),params.secset);
else
    secstr = '';
end

Nconj = sum(ndx);

logstr = [current_timestring() ...
    ' Conjunctions involving the primary' secstr ' set: ' ...
    num2str(Nconj) ' of ' num2str(Ntot)];
log_string(params.logfid,logstr,params.logging,params.displaying);

% Restrict output DB to the specified primaries
out.ConjDist.DB = out.ConjDist.DB(ndx,:);

% Return for no conjunctions involving the primary set
if (Nconj == 0)
    if params.logging; fclose(params.logfid); end
    return;
end

%% Associate updates into "events" (i.e., conjunction update sequences)
% which all have the same primary and secondary and very nearly the same
% TCA

TCAdel = 0.5*(params.TCAtol_s/86400);

out.ConjDist.event = NaN(Nconj,1);
Nevent = 0;

for n=1:Nconj
    
    % Only attempt to associate previously unassociated conjunctions
    
    if isnan(out.ConjDist.event(n))
        
        % Secondary
        
        p = out.ConjDist.DB(n,1);
        s = out.ConjDist.DB(n,2);
        
        % TCA (days since January 0, 0000)
        
        TCA = out.ConjDist.DB(n,166);
        
        TCA1 = TCA-TCAdel;
        TCA2 = TCA+TCAdel;
        
        % Associate conjunctions
        
        ndx = isnan(out.ConjDist.event)       & ...
              (p == out.ConjDist.DB(:,1))     & ...
              (s == out.ConjDist.DB(:,2))     & ...
              (out.ConjDist.DB(:,166) >= TCA1 & ...
              (out.ConjDist.DB(:,166) <= TCA2));
        Nndx = sum(ndx);
        
        if (Nndx == 0)
            error('Failure to find any updates');
        end

        Nevent = Nevent + 1;
        out.ConjDist.event(ndx) = Nevent;
        
        % Check for disagreements with OCMDB file associations
        
        if params.association_check
            origevent = out.ConjDist.DB(n,217);
            mdx = (out.ConjDist.DB(:,217) == origevent);
            mdx_ne_ndx = ~isequal(ndx,mdx);
            if mdx_ne_ndx
                wrnstr = ['Association disagreement: n = ' num2str(n) ...
                    ' p = ' num2str(p) ' s = ' num2str(s) ...
                    ' TCA = ' datestr(TCA,31)];
                warning(wrnstr);
                logstr = ['WARNING: ' wrnstr];
                log_string(params.logfid,logstr,params.logging,params.displaying);
                for npass=1:2
                    if npass == 1
                        idx = find(ndx);
                        Nidx = numel(idx);
                        logstr = ['  ConjDist associations:'];
                    else
                        idx = find(mdx);
                        Nidx = numel(idx);
                        logstr = ['  OCMDB associations:'];
                    end
                    log_string(params.logfid,logstr,params.logging,params.displaying);
                    mdT = median(out.ConjDist.DB(idx,166));
                    mdX = median(out.ConjDist.DB(idx,172));
                    mdY = median(out.ConjDist.DB(idx,173));
                    mdZ = median(out.ConjDist.DB(idx,174));
                    for nn=1:Nidx
                        ii = idx(nn);
                        logstr = ['   ' datestr(out.ConjDist.DB(ii,166),'yyyy-mm-dd HH:MM:SS.FFF') ...
                            ' Dt=' num2str(86400*(mdT-out.ConjDist.DB(ii,166))) 's ' ...
                            ' Dx=' num2str(mdX-out.ConjDist.DB(ii,172)) 'm ' ...
                            ' Dy=' num2str(mdY-out.ConjDist.DB(ii,173)) 'm ' ...
                            ' Dz=' num2str(mdZ-out.ConjDist.DB(ii,174)) 'm '];
                        log_string(params.logfid,logstr,params.logging,params.displaying);
                    end
                end
            end
        end
            
    end
    
end

if any(isnan(out.ConjDist.event))
    error('Failure to associate all conjunction updates');
end

out.ConjDist.Nevent = Nevent;

logstr = [current_timestring() ...
    ' Associated unique events involving the primary set: ' ...
    num2str(Nevent)];
log_string(params.logfid,logstr,params.logging,params.displaying);

% Return for no conjunctions involving the primary set

if (Nevent == 0)
    if params.logging; fclose(params.logfid); end
    return;
end

% Calculate the median TCA for each event

out.ConjDist.TCAmed = NaN(Nevent,1);

for ne=1:Nevent
    ndx = (ne == out.ConjDist.event);
    out.ConjDist.TCAmed(ne) = median(out.ConjDist.DB(ndx,166));
end

minTCAmed = min(out.ConjDist.TCAmed);
maxTCAmed = max(out.ConjDist.TCAmed);

minTCA = round_timestring(jd_to_timestring(JD0+minTCAmed),'m');
maxTCA = round_timestring(jd_to_timestring(JD0+maxTCAmed),'m');
delTCA_yr = (maxTCAmed-minTCAmed)/365.25;

logstr = [current_timestring() ...
    ' Time spanned by median event TCA values: '];
log_string(params.logfid,logstr,params.logging,params.displaying);
logstr = [pad20 minTCA ' to ' maxTCA ' (' ...
    num2str(delTCA_yr,'%0.2f') ' yr)'];
log_string(params.logfid,logstr,params.logging,params.displaying);

% Eliminate events with median TCAs outside the TCA limits

out.ConjDist.Nconj0 = Nconj;
out.ConjDist.Nevent0 = out.ConjDist.Nevent;

if impose_TCAlimits
    
    logstr = [current_timestring() ...
        ' Restricting to median event TCA limits: '];
    log_string(params.logfid,logstr,params.logging,params.displaying);
    if strcmpi(class(params.TCAlimits_UT),'double')
        logstr = ['  ' num2str(params.TCAlimits_UT(1)) ' days after beginning'...
                ' to ' num2str(params.TCAlimits_UT(2)) ' days before ending'];
    else
        logstr = ['  ' params.TCAlimits_UT{1} ' to ' params.TCAlimits_UT{2}];
    end
    log_string(params.logfid,logstr,params.logging,params.displaying);

    % Find events within the TCA limits

    inlimits = find( (out.ConjDist.TCAlimits(1) <= out.ConjDist.TCAmed) & ...
                     (out.ConjDist.TCAlimits(2) >= out.ConjDist.TCAmed) );
    Ninlimits = numel(inlimits);
    
    if (Ninlimits < Nevent)
        
        % Some events are not within TCA limits and need to be excluded
        
        % Find the index values for all raw conjunctions within the TCA
        % limits
        
        for ni=1:Ninlimits
            ne = inlimits(ni);
            if (ni == 1)
                ndx = (ne == out.ConjDist.event);
            else
                ndx = ndx | (ne == out.ConjDist.event);
            end
        end
        
        % Update the number of raw conjunctions, and restrict the output DB
        % to those within the TCA limits
        
        Nconj = sum(ndx);
        out.ConjDist.DB = out.ConjDist.DB(ndx,:);
        
        % Update the events and their median TCA values

        resevent = out.ConjDist.event(ndx);
        unqevent = unique(resevent);
        Nunqevent = numel(unqevent);
        newevent = NaN(Nconj,1);
        for nue=1:Nunqevent
            idx = (unqevent(nue) == resevent);
            newevent(idx) = nue;
        end
        
        Nevent = Ninlimits;
        out.ConjDist.Nevent = Nevent;
        out.ConjDist.event = newevent;        
        out.ConjDist.TCAmed = out.ConjDist.TCAmed(inlimits);

    end
    
    logstr = [current_timestring() ...
        ' Associated conjunction events within TCA limits: ' ...
        num2str(Nevent)];
    log_string(params.logfid,logstr,params.logging,params.displaying);
    
    minTCAmed = min(out.ConjDist.TCAmed);
    maxTCAmed = max(out.ConjDist.TCAmed);

    minTCA = round_timestring(jd_to_timestring(JD0+minTCAmed),'m');
    maxTCA = round_timestring(jd_to_timestring(JD0+maxTCAmed),'m');
    delTCA_yr = (maxTCAmed-minTCAmed)/365.25;    

    logstr = [current_timestring() ...
        ' Time spanned by median event TCA values within limits: '];
    log_string(params.logfid,logstr,params.logging,params.displaying);
    logstr = ['  ' minTCA ' to ' maxTCA ' (' ...
        num2str(delTCA_yr,'%0.2f') ' yr)'];        
    log_string(params.logfid,logstr,params.logging,params.displaying);
    
end

out.ConjDist.Nconj = Nconj;

%% Apogee, perigee and inclination extracted from DB array

ii = 51;
mdapo = median(out.ConjDist.DB(:,ii));
mnapo = mean(out.ConjDist.DB(:,ii));
sgapo = std(out.ConjDist.DB(:,ii));

ii = 52;
mdper = median(out.ConjDist.DB(:,ii));
mnper = mean(out.ConjDist.DB(:,ii));
sgper = std(out.ConjDist.DB(:,ii));

ii = 53;
mdinc = median(out.ConjDist.DB(:,ii));
mninc = mean(out.ConjDist.DB(:,ii));
sginc = std(out.ConjDist.DB(:,ii));

logstr = [current_timestring() ' Primary set perigee median, mean, std: ' ...
    num2str(mdper) ' ' num2str(mnper) ' ' num2str(sgper)];
log_string(params.logfid,logstr,params.logging,params.displaying);

logstr = [current_timestring() ' Primary set apogee median, mean, std: ' ...
    num2str(mdapo) ' ' num2str(mnapo) ' ' num2str(sgapo)];
log_string(params.logfid,logstr,params.logging,params.displaying);

logstr = [current_timestring() ' Primary set inclination median, mean, std: ' ...
    num2str(mdinc) ' ' num2str(mninc) ' ' num2str(sginc)];
log_string(params.logfid,logstr,params.logging,params.displaying);

%% Find all of the unique secondaries and recalculate Pc values if required

% Unique secondaries
sec = unique(out.ConjDist.DB(:,2));
out.ConjDist.Nsec = numel(sec);

logstr = [current_timestring() ...
    ' Unique secondaries in associated conjunction events: ' ...
    num2str(out.ConjDist.Nsec)];
log_string(params.logfid,logstr,params.logging,params.displaying);

% Allocate array for new Pc values
Pcnew = NaN(Nconj,1);

% Use Pc interpolation table, if available, otherwise calculate values
if ~using_PcTable

    % Use OCMDB states and covariances to calculate new Pc values

    logstr = [current_timestring() ...
        ' Calculating ' num2str(Nconj) ' Pc values' ...
        ' with HBRnew = ' num2str(params.mission_HBR_meters) 'm)'];
    log_string(params.logfid,logstr,params.logging,params.displaying);

    % Extract stats in meter units
    DBX1 = out.ConjDist.DB(:,172:177)*1e3;
    DBC1 = out.ConjDist.DB(:,73:93);
    DBX2 = out.ConjDist.DB(:,178:183)*1e3;
    DBC2 = out.ConjDist.DB(:,133:153);
    HBRmeters = params.mission_HBR_meters;

    hbar = parfor_progressbar(Nconj,'','Name','New HBR Pc Calculation');

    for nc=1:Nconj % parfor enabled
        % Primary object ECI state
        X1 = DBX1(nc,:);
        % Primary object covariance 
        TriEls = DBC1(nc,:);
        C1UVW = TriEls(CovBuild6);
        % Convert UVW (same as the RIC) to ECI
        C1ECI = RIC2ECI(C1UVW,X1(1:3),X1(4:6));
        C1ECI = cov_make_symmetric(C1ECI);
        % Secondary object ECI state
        X2 = DBX2(nc,:);
        % Secondary object covariance 
        TriEls = DBC2(nc,:);
        C2UVW = TriEls(CovBuild6);
        % Convert UVW (same as the RIC) to ECI
        C2ECI = RIC2ECI(C2UVW,X2(1:3),X2(4:6));
        C2ECI = cov_make_symmetric(C2ECI);
        % Calculate Pc
        Pcnew(nc) = PcMultiStep( ...
            X1(1:3),X1(4:6),C1ECI, ...
            X2(1:3),X2(4:6),C2ECI, ...
            HBRmeters,PcMultiStep_params);

        hbar.iterate(1);

    end

    clear DBX1 DBC1 DBX2 DBC2;
    close(hbar);        

else

    % Use Pc interpolation table to calculate new Pc values

    logstr = [current_timestring() ...
        ' Interpolating ' num2str(Nconj) ' new Pc ' ...
        'values because of potential HBR mismatch (HBRnew=' ...
        num2str(params.mission_HBR_meters) 'm)'];
    log_string(params.logfid,logstr,params.logging,params.displaying);

    % Ensure HBR is within the Pc table range
    if params.PcTable.HBRmin_meters > params.mission_HBR_meters || ...
       params.PcTable.HBRmax_meters < params.mission_HBR_meters
       error('Specified mission HBR must be within PcTable range');
    else
        % First determine if HBR is on the PcTable grid
        [~,nHBR0] = min(abs(params.PcTable.HBR-params.mission_HBR_meters));
        if params.PcTable.HBR(nHBR0) == params.mission_HBR_meters
            HBR_on_grid = true;
        else
            % Define points bracketing HBR in table        
            HBR_on_grid = false;
            log10HBRnew = log10(params.mission_HBR_meters);
            nHBR1 = nHBR0-3;
            if nHBR1 < 1
                nHBR1 = 1; nHBR2 = nHBR1+6;
            else
                nHBR2 = nHBR0+3;
                if nHBR2 > params.PcTable.HBRnum
                    nHBR2 = params.PcTable.HBRnum; nHBR1 = nHBR2-6;
                end
            end
        end
    end

    % Get or interpolate the Pc value for the new HBR from the Pc table
    Ncol0 = params.PcTable.OriginalDBColumns;
    if HBR_on_grid
        % Get Pc values from table previously stored in DB array
        Pcnew(:) = out.ConjDist.DB(:,Ncol0+nHBR0);
    else
        % Interpolate Pc values from table previously stored in DB
        nHBR = nHBR1:nHBR2;
        iHBR = Ncol0+nHBR;
        for nc=1:Nconj
            % Linear interpolation if any of the local grid Pc values
            % are zero, otherwise logarithmic interpolation
            if any(out.ConjDist.DB(nc,iHBR) == 0)
                Pcnew(nc) = interp1(params.PcTable.HBR(nHBR), ...
                                    out.ConjDist.DB(nc,iHBR), ...
                                    params.mission_HBR_meters);
            else
                Pcnew(nc) = 10^interp1(params.PcTable.log10HBR(nHBR), ...
                                    log10(out.ConjDist.DB(nc,iHBR)), ...
                                    log10HBRnew);
            end
        end
    end

end

% Insert the new HBR and Pc into the DB array
out.ConjDist.DB(:,157) = params.mission_HBR_meters;
out.ConjDist.DB(:,156) = Pcnew;
clear Pcnew;

%% Create a list of events

logstr = [current_timestring() ...
    ' Defining output list of associated conjunction events'];
log_string(params.logfid,logstr,params.logging,params.displaying);

% Open and initialize CSV output file listing all associated events

if params.make_event_file
    
    csv_file_root = ['Events_' params.prisetstr '_Ne' num2str(Nevent)];
    if params.timestamp
        csv_file_root = [csv_file_root '_' params.timetag];
    end
    csv_file_name = fullfile(params.outputdir,[csv_file_root '.csv']);
    csvfid = fopen(csv_file_name,'wt');
    csvfmt = '%0.17g';

    outstr = 'Primary';
    outstr = [outstr ', ' 'Secondary'];
    outstr = [outstr ', ' 'Nupdate'];
    outstr = [outstr ', ' 'medTCA'];
    outstr = [outstr ', ' 'maxdelTCA'];
    if ~using_PcTable
        outstr = [outstr ', ' 'NCAinSV'];
        outstr = [outstr ', ' 'NLTinSV'];
    end
    outstr = [outstr ', ' 'minPc'];
    outstr = [outstr ', ' 'maxPc'];
    outstr = [outstr ', ' 'medPc'];
    outstr = [outstr ', ' 'fudPc'];
    outstr = [outstr ', ' 'ludPc'];
    outstr = [outstr ', ' 'ludAge'];
    outstr = [outstr ', ' 'ludUse'];
    outstr = [outstr ', ' 'mudPc'];
    outstr = [outstr ', ' 'mudAge'];
    outstr = [outstr ', ' 'mudUse'];
    outstr = [outstr ', ' 'minHBR'];
    outstr = [outstr ', ' 'maxHBR'];

    fprintf(csvfid,'%s\n',outstr);
    
end

% Allocate arrays for event processing
out.ConjDist.Pri    = NaN(Nevent,1);
out.ConjDist.Sec    = NaN(Nevent,1);
out.ConjDist.TCAfud = NaN(Nevent,1);
out.ConjDist.TCAlud = NaN(Nevent,1);
out.ConjDist.TCAmud = NaN(Nevent,1);
out.ConjDist.TCAmin = NaN(Nevent,1);
out.ConjDist.TCAmax = NaN(Nevent,1);
out.ConjDist.Pcmin  = NaN(Nevent,1);
out.ConjDist.Pcmax  = NaN(Nevent,1);
out.ConjDist.Pcmed  = NaN(Nevent,1);
out.ConjDist.Pcfud  = NaN(Nevent,1);
out.ConjDist.Pclud  = NaN(Nevent,1);
out.ConjDist.Pcmud  = NaN(Nevent,1);
out.ConjDist.AGElud = NaN(Nevent,1);
out.ConjDist.USElud = NaN(Nevent,1);
out.ConjDist.AGEmud = NaN(Nevent,1);
out.ConjDist.USEmud = NaN(Nevent,1);
out.ConjDist.HBRmin = NaN(Nevent,1);
out.ConjDist.HBRmax = NaN(Nevent,1);
if ~using_PcTable
    out.ConjDist.NCAin  = zeros(Nevent,1);
    out.ConjDist.NLTin  = zeros(Nevent,1);
end

% Parameters usage of CArate_in_screen function

clear CApr;
CApr.screen_type = params.SVRICtype_alt;

% Perform event processing

N_tinyRCSorSpeed = 0;
N_noncatastrophic = 0;
M_noncatastrophic = [];

for ne=1:Nevent
    
    % Indices for this event's update sequence
    
    ndx0 = (ne == out.ConjDist.event);
    ndx  = find(ndx0);
    Nndx = numel(ndx);
    
    % Primary and secondary for this event
    
    p = out.ConjDist.DB(ndx(1),1); out.ConjDist.Pri(ne) = p;
    s = out.ConjDist.DB(ndx(1),2); out.ConjDist.Sec(ne) = s;
    
    % Index array for this primary
    
    pdx = (params.priset == p);
    
    % Parameters array for usage of CArate_in_screen for screening volume
    % incursions

    CApr.screen_dimensions = params.SVRICdims_m(pdx,:);
    CApr.screen_max_r = max(CApr.screen_dimensions);
    CApr.screen_max_r2 = CApr.screen_max_r^2;
    CApr.screen_ratios = CApr.screen_dimensions / ...
                         CApr.screen_max_r;
    CApr.Cnorm = diag(CApr.screen_ratios.^(-2));
    
    if ~strcmpi(params.SVRICtype,'box')
        aaa = params.SVRICdims_m(pdx,1);
        bbb = params.SVRICdims_m(pdx,2);
        ccc = params.SVRICdims_m(pdx,3);
        aaa2 = aaa^2;
        bbb2 = bbb^2;
        ccc2 = ccc^2;
    end
    
    % TCA for this events update sequence
    
    TCAndx = out.ConjDist.DB(ndx,166);
    TCAmedcheck = median(TCAndx);
    if (out.ConjDist.TCAmed(ne) ~= TCAmedcheck)
        error('Median TCA mismatch');
    end
    
    out.ConjDist.TCAmin(ne) = min(TCAndx);
    out.ConjDist.TCAmax(ne) = max(TCAndx);
    
    TCRndx = out.ConjDist.DB(ndx,165);
    HBRndx = out.ConjDist.DB(ndx,157);
    if any(HBRndx ~= params.mission_HBR_meters)
        error('HBR mismatch found');
    end
        
    % Calculate Pc values
    Pcndx = out.ConjDist.DB(ndx,156);
    
    % Clip small Pc values
    Pcndx(Pcndx <= 1e-300) = 0;
    
    % Min/max HBR values
    out.ConjDist.HBRmin(ne) = min(HBRndx);
    out.ConjDist.HBRmax(ne) = max(HBRndx);
    
    % First Pc update in sequence
    [~,imin] = min(TCRndx);
    out.ConjDist.TCAfud(ne) = TCAndx(imin);
    out.ConjDist.Pcfud(ne)  = Pcndx(imin);
    
    % Find conjunctions in Pc update sequence within the commit and
    % consider times
    TCAmTCR = TCAndx-TCRndx;
    InCommitConsider = (TCAmTCR <= max(params.CommitConsiderDays)) & ...
                       (TCAmTCR >= min(params.CommitConsiderDays));

    if any(InCommitConsider)
        
        % Find updates found within specified limits
        InCommitConsider = find(InCommitConsider);
        
        % Indicate this event's last-update and max-update Pc values can be
        % used because they are within the consider/commit time limits
        out.ConjDist.USElud(ne) = 1;
        out.ConjDist.USEmud(ne) = 1;
        
        % Index of last-update in this event's conjunction list
        [out.ConjDist.AGElud(ne),ilud] = min(TCAmTCR(InCommitConsider));
        ilud = InCommitConsider(ilud);
        
        % Index of max Pc value in this events conjunction list
        [~,imud] = max(Pcndx(InCommitConsider));
        imud = InCommitConsider(imud);
        out.ConjDist.AGEmud(ne) = TCAmTCR(imud);

        % Check for catastrophic status of last-update conjunction, if
        % required
        if params.exclude_noncatastrophic
            if PcLUDvsPcMaxOption == 1
                nc = ndx(ilud); % Index of last-update in DB conjunction list
            elseif PcLUDvsPcMaxOption == 2
                nc = ndx(imud); % Index of max Pc in DB conjunction list
            else
                error('Invalid PcLUDvsPcMaxOption for exclude_noncatastrophic = true');
            end
            tinyRCS = 1e-12;
            SecRCS = out.ConjDist.DB(nc,114);
            RelSpeed = out.ConjDist.DB(nc,20);
            if (SecRCS < tinyRCS) || (RelSpeed <= 0)
                % Zero or tiny secondary RCS values cannot be used to
                % evaluate chance of non-catastophic collision
                N_tinyRCSorSpeed = N_tinyRCSorSpeed+1;
            else
                % Estimate catastophic events using collision consequence
                % method
                [Catastrophic,~,massVec] = ...
                    CollisionConsequenceOCMDBReader(out.ConjDist.DB(nc,:), ...
                        params.mission_on_orbit_mass_kg);
                ChanceNonCatastrophic = 1-sum(Catastrophic)/length(Catastrophic);
                if ChanceNonCatastrophic > params.confidence_noncatastrophic
                    % Do not use high-confidence noncatastrophic collisions
                    if PcLUDvsPcMaxOption == 1
                        out.ConjDist.USElud(ne) = 0;
                    elseif PcLUDvsPcMaxOption == 2
                        out.ConjDist.USEmud(ne) = 0;
                    end
                    N_noncatastrophic = N_noncatastrophic+1;
                    M_noncatastrophic = cat(1,M_noncatastrophic,median(massVec));
                end
            end
        end
        
    else
        
        % No updates found within specified limits
        % This events's last-update Pc can't be used
        out.ConjDist.USElud(ne) = 0; 
        [out.ConjDist.AGElud(ne),ilud] = min(TCAmTCR);
        % This events's max Pc can't be used
        out.ConjDist.USEmud(ne) = 0;
        [~,imud] = max(Pcndx);
        out.ConjDist.AGEmud(ne) = TCAmTCR(imud);
        
    end

    % Last update and max-Pc update Pc values and TCAs
    out.ConjDist.Pclud(ne)  = Pcndx(ilud);
    out.ConjDist.TCAlud(ne) = TCAndx(ilud);
    out.ConjDist.Pcmud(ne)  = Pcndx(imud);
    out.ConjDist.TCAmud(ne) = TCAndx(imud);
    
    % Median TCA string
    TCAmedstr = round_timestring(jd_to_timestring(JD0+out.ConjDist.TCAmed(ne)),'s');
    
    dTCAmax = max(abs([out.ConjDist.TCAmin(ne)-out.ConjDist.TCAmed(ne) ...
                       out.ConjDist.TCAmax(ne)-out.ConjDist.TCAmed(ne)]));
    dTCAmaxstr = num2str(86400*dTCAmax,'%+0.3f');

    out.ConjDist.Pcmin(ne) = min(Pcndx);
    out.ConjDist.Pcmax(ne) = max(Pcndx);
    out.ConjDist.Pcmed(ne) = median(Pcndx);
    
    % Calculate RIC screening volume penetrations for all conjunction updates
    
    if ~using_PcTable
        
        % Use state info in DB to determine number of close-approach points
        % and linear-trajectories that are within the screening volume
        for nn=1:Nndx

            % Index for this conjunction update

            n = ndx(nn);

            % Get the secondary's state in the primary's RIC frame at TCA.

            % Calculate the secondary's position in the primary's RIC frame

            pos1 = out.ConjDist.DB(n,172:174)*1e3; % ECI primary pos & vel
            vel1 = out.ConjDist.DB(n,175:177)*1e3;

            pos2 = out.ConjDist.DB(n,178:180)*1e3; % TDR secondary pos & vel
            vel2 = out.ConjDist.DB(n,181:183)*1e3;

            h = cross(pos1,vel1);
            rhat = pos1 / norm(pos1);
            chat = h / norm(h);
            ihat = cross(chat,rhat);
            ECItoRIC = [rhat; ihat; chat];

            % Cross check for low- vs high-precision

            rs = ECItoRIC * (pos2-pos1)';
            vs = ECItoRIC * (vel2-vel1)';

            % Determine if the CA point and/or the linearized trajectory (LT)
            % penetrates the RIC screening volume

            [CA_in_SV, LT_in_SV] = CArate_in_screen( ...
                CApr.screen_max_r, CApr.screen_max_r2, true, ...
                rs, vs, CApr);

            % Cross check for ellipsoidal SV

            if ~strcmpi(params.SVRICtype,'box')

                % Ellipsoid cross-check code 1

                qqq1 = rs;
                sca1 = qqq1(1)/aaa; sca2 = qqq1(2)/bbb; sca3 = qqq1(3)/ccc;
                fca = sca1^2 + sca2^2 + sca3^2 - 1;
                CAin = fca <= 0;

                ttop = ( rs(1)*vs(1)/aaa2 + ...
                         rs(2)*vs(2)/bbb2 + ...
                         rs(3)*vs(3)/ccc2 );

                tbot = ( vs(1)*vs(1)/aaa2 + ...
                         vs(2)*vs(2)/bbb2 + ...
                         vs(3)*vs(3)/ccc2 );

                ttt0 = -ttop/tbot;
                qqq0 = rs + ttt0*vs;
                slt1 = qqq0(1)/aaa; slt2 = qqq0(2)/bbb; slt3 = qqq0(3)/ccc;
                flt = slt1^2 + slt2^2 + slt3^2 - 1;
                LTin = flt <= 0;
                if (CA_in_SV ~= CAin) || (LT_in_SV ~= LTin)
                    % keyboard;
                    error('Screening ellipsoid cross check 1 failed');
                end

                % Ellipsoid cross-check code 2

                sigRIC = [aaa; bbb; ccc];
                rr = rs./sigRIC;
                vv = vs./sigRIC;
                AA = sum(vv.^2);
                BB = 2*sum(rr.*vv);
                CC = sum(rr.^2)-1;
                ca_in_sv = CC < 0;
                lt_in_sv = BB^2-4*AA*CC > 0;
                if (ca_in_sv ~= CA_in_SV) || (lt_in_sv ~= LT_in_SV)
                    % keyboard;
                    error('Screening ellipsoid cross check 1 failed');
                end

            end

            out.ConjDist.NCAin(ne) = out.ConjDist.NCAin(ne) + CA_in_SV;
            out.ConjDist.NLTin(ne) = out.ConjDist.NLTin(ne) + LT_in_SV;

        end

    end
    
    % Output results
    
    if params.make_event_file    
        
        outstr = num2str(p,'%05i');
        outstr = [outstr ', ' num2str(s,'%05i')];
        outstr = [outstr ', ' num2str(Nndx)];
        outstr = [outstr ', ' TCAmedstr];
        outstr = [outstr ', ' dTCAmaxstr];

        if ~using_PcTable
            outstr = [outstr ', ' num2str(out.ConjDist.NCAin(ne))];
            outstr = [outstr ', ' num2str(out.ConjDist.NLTin(ne))];
        end
        outstr = [outstr ', ' num2str(out.ConjDist.Pcmin(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.Pcmax(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.Pcmed(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.Pcfud(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.Pclud(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.AGElud(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.USElud(ne))];
        outstr = [outstr ', ' num2str(out.ConjDist.Pcmud(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.AGEmud(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.USEmud(ne))];
        outstr = [outstr ', ' num2str(out.ConjDist.HBRmin(ne),csvfmt)];
        outstr = [outstr ', ' num2str(out.ConjDist.HBRmax(ne),csvfmt)];
        
        fprintf(csvfid,'%s\n',outstr);
        
    end
    
    % Make update plot for this event if required
    
    make_event_update_plot = false;
    if params.make_update_plots
        if PcLUDvsPcMaxOption == 1
            % For latest update option, ensure the LUD is within the usable
            % age range, and the LUD Pc value is sufficiently large
            make_event_update_plot = ...
                (out.ConjDist.USElud(ne) == 1) & ...
                (out.ConjDist.Pcmax(ne) >= params.Pc_for_update_plots);
        elseif PcLUDvsPcMaxOption == 2
            % For max-Pc option, ensure the MUD is within the usable
            % age range, and the MUD Pc value is sufficiently large
            make_event_update_plot = ...
                (out.ConjDist.USEmud(ne) == 1) & ...
                (out.ConjDist.Pcmax(ne) >= params.Pc_for_update_plots);
        else
            % Ensure that there are enough updates and that the max Pc
            % value is sufficiently large
            make_event_update_plot = ...
                (Nndx >= params.Nu_for_update_plots) & ...
                (out.ConjDist.Pcmax(ne) >= params.Pc_for_update_plots);
        end
    end
    
    if make_event_update_plot

        % ID strings for earliest and latest updates

        [~,jmin] = min(TCRndx);
        nc = ndx(jmin);
        id1 = sprintf(['%05i%s%05i%s%04i%02i%02i%s%02i' ...
            '%02i%02i%s%04i%02i%02i%s%02i%02i%02i'], ...
            out.ConjDist.DB(nc,1),'_conj_',out.ConjDist.DB(nc,2),'_', ...
            out.ConjDist.DB(nc,12:14),'_',out.ConjDist.DB(nc,15:17),'_', ...
            out.ConjDist.DB(nc,3:5),'_',out.ConjDist.DB(nc,6:8));

        [~,jmax] = max(TCRndx);
        nc = ndx(jmax);
        id2 = sprintf(['%05i%s%05i%s%04i%02i%02i%s%02i' ...
            '%02i%02i%s%04i%02i%02i%s%02i%02i%02i'], ...
            out.ConjDist.DB(nc,1),'_conj_',out.ConjDist.DB(nc,2),'_', ...
            out.ConjDist.DB(nc,12:14),'_',out.ConjDist.DB(nc,15:17),'_', ...
            out.ConjDist.DB(nc,3:5),'_',out.ConjDist.DB(nc,6:8));
        
        % Construct title
        
        nti = 0; clear titl;
        nti=nti+1; titl{nti} = 'Collision Probability Updates';
        nti=nti+1; titl{nti} = strrep(id1,'_','\_');
        nti=nti+1; titl{nti} = strrep(id2,'_','\_');
        nti=nti+1; titl{nti} = ['Pc values estimated using HBR = ' ...
                                num2str(params.mission_HBR_meters) ' m'];
        if PcLUDvsPcMaxOption == 1
            nti=nti+1;
            titl{nti} = ['Last update ' num2str(min(params.CommitConsiderDays)) ...
                ' to ' num2str(max(params.CommitConsiderDays)) ...
                ' days before TCA: Pc = ' smart_exp_format(out.ConjDist.Pclud(ne),3)];
        elseif PcLUDvsPcMaxOption == 2
            nti=nti+1;
            titl{nti} = ['Max-Pc update ' num2str(min(params.CommitConsiderDays)) ...
                ' to ' num2str(max(params.CommitConsiderDays)) ...
                ' days before TCA: Pc = ' smart_exp_format(out.ConjDist.Pcmud(ne),3)];
        end
        nti=nti+1; titl{nti} = ' ';
        
        % Make update plot
        
        fighandle1 = ConjDist_update_plot(titl, ...
            out.ConjDist.TCAlud(ne),TCRndx,Pcndx,params,fighandle1);
        
        plt_file_root = [id2 '_Nupdate' num2str(Nndx)];        
        plt_file_name = fullfile(params.outputdir,[plt_file_root params.plot_format]);
        saveas(fighandle1,plt_file_name);

    end
    
end

if params.exclude_noncatastrophic && (N_noncatastrophic > 0)
    logstr = [current_timestring() ...
        ' (SecRCS < ' smart_exp_format(tinyRCS,3) ...
        ') or (RelSpeed <= 0) events counted: ' ...
        num2str(N_tinyRCSorSpeed) ' out of ' num2str(Nevent)];
    log_string(params.logfid,logstr,params.logging,params.displaying);
    logstr = ['  ' num2str(params.confidence_noncatastrophic*100) '%-likely' ...
        ' noncatastrophic events counted: ' num2str(N_noncatastrophic) ...
        ' out of ' num2str(Nevent)];
    log_string(params.logfid,logstr,params.logging,params.displaying);
    if N_noncatastrophic > 0
        logstr = ['  Min/med/max sec. mass for likely non-cat. events = ' ...
            num2str(min(M_noncatastrophic)) '/' ...
            num2str(median(M_noncatastrophic)) '/' ...
            num2str(max(M_noncatastrophic)) ' kg'];
        log_string(params.logfid,logstr,params.logging,params.displaying);
    end
end

% Close the figure and the csv file

if ishandle(fighandle1); delete(fighandle1); end

if params.make_event_file; fclose(csvfid); end

% Make repeating event file (turned off for EventRate version)

if params.make_repeating_event_file
    % Generate the repeating events
    [repev] = find_repeating_events(csv_file_name,out.ConjDist.DB,params);
end

% Issue warning if any HBR changes were detected

HBRminall = min(out.ConjDist.HBRmin);
HBRmaxall = max(out.ConjDist.HBRmax);

if (HBRminall ~= HBRmaxall)
    wrnstr = ['HBR changes occurred during time span:' ...
        ' minHBR = ' num2str(HBRminall) ...
        ' maxHBR = ' num2str(HBRmaxall)];
    logstr = [current_timestring() ...
        ' WARNING: ' wrnstr];
    log_string(params.logfid,logstr,params.logging,params.displaying);
    warning(wrnstr);
end

%% Time span covered by the associated conjunction events

out.ConjDist.TCAmin = minTCAmed;
out.ConjDist.TCAmax = maxTCAmed;
out.ConjDist.TCAdelt = maxTCAmed-minTCAmed;

logstr = [current_timestring() ....
    ' Time spanned by the associated event median TCA values: ' ...
    num2str(out.ConjDist.TCAdelt) ' days'];
log_string(params.logfid,logstr,params.logging,params.displaying);

% Generate the output TCA bins

if params.make_time_bins
    % Bin width
    Wbin = params.TCAbin_days;
    % Calculate the number of bins. Use floor here to prevent the
    % possiblity of an underfilled last bin
    out.ConjDist.bins.Nbin = floor((out.ConjDist.TCAlimits(2)-out.ConjDist.TCAlimits(1))/Wbin);
    out.ConjDist.bins.Tbinlo = out.ConjDist.TCAlimits(1)+params.TCAbin_days*(0:out.ConjDist.bins.Nbin-1);
    out.ConjDist.bins.Tbinhi = out.ConjDist.bins.Tbinlo+Wbin;
    out.ConjDist.bins.Tbinhi(out.ConjDist.bins.Nbin) = ...
        max(out.ConjDist.bins.Tbinhi(out.ConjDist.bins.Nbin),out.ConjDist.TCAlimits(2));
    out.ConjDist.bins.Tbin = (out.ConjDist.bins.Tbinhi+out.ConjDist.bins.Tbinlo)/2;
end

%% Report and plot all registered events and yellow/red events, as required
%
% The can comprise three sets:
%   1) the entire set of associated conjunction events
%   2) the subset representing primary SVI incursions
%   3) the complementary subset not representing primary SVI incursions
%
% The EventRate version of ConjDist only uses set 1, and the other two are
% just turned off, rather than eliminated

out.ConjDist.Nall = NaN(3,1);
out.ConjDist.Nyel = NaN(3,1);
out.ConjDist.Nred = NaN(3,1);

out.ConjDist.rates.all.md = NaN(3,1);
out.ConjDist.rates.all.lo = NaN(3,1);
out.ConjDist.rates.all.hi = NaN(3,1);

out.ConjDist.rates.yel.md = NaN(3,1);
out.ConjDist.rates.yel.lo = NaN(3,1);
out.ConjDist.rates.yel.hi = NaN(3,1);

out.ConjDist.rates.red.md = NaN(3,1);
out.ConjDist.rates.red.lo = NaN(3,1);
out.ConjDist.rates.red.hi = NaN(3,1);

out.ConjDist.fracs.yel.md = NaN(3,1);
out.ConjDist.fracs.yel.lo = NaN(3,1);
out.ConjDist.fracs.yel.hi = NaN(3,1);

out.ConjDist.fracs.red.md = NaN(3,1);
out.ConjDist.fracs.red.lo = NaN(3,1);
out.ConjDist.fracs.red.hi = NaN(3,1);

if params.make_time_bins
    
    out.ConjDist.rates.binall.md = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binall.lo = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binall.hi = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binyel.md = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binyel.lo = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binyel.hi = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binred.md = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binred.lo = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.rates.binred.hi = NaN(3,out.ConjDist.bins.Nbin);
    
    out.ConjDist.fracs.binyel.md = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.fracs.binyel.lo = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.fracs.binyel.hi = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.fracs.binred.md = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.fracs.binred.lo = NaN(3,out.ConjDist.bins.Nbin);
    out.ConjDist.fracs.binred.hi = NaN(3,out.ConjDist.bins.Nbin);
    
    out.ConjDist.rates.bavgall.md = NaN(3,1);
    out.ConjDist.rates.bavgall.lo = NaN(3,1);
    out.ConjDist.rates.bavgall.hi = NaN(3,1);
    out.ConjDist.rates.bavgyel.md = NaN(3,1);
    out.ConjDist.rates.bavgyel.lo = NaN(3,1);
    out.ConjDist.rates.bavgyel.hi = NaN(3,1);
    out.ConjDist.rates.bavgred.md = NaN(3,1);
    out.ConjDist.rates.bavgred.lo = NaN(3,1);
    out.ConjDist.rates.bavgred.hi = NaN(3,1);

    out.ConjDist.fracs.bavgyel.md = NaN(3,1);
    out.ConjDist.fracs.bavgyel.lo = NaN(3,1);
    out.ConjDist.fracs.bavgyel.hi = NaN(3,1);
    out.ConjDist.fracs.bavgred.md = NaN(3,1);
    out.ConjDist.fracs.bavgred.lo = NaN(3,1);
    out.ConjDist.fracs.bavgred.hi = NaN(3,1);
    
end

% If last-update option, then only process full set - not those that are
% screening volume incursions (SVI), or non-SVIs
if PcLUDvsPcMaxOption > 0
    nset2 = 1;
else
    if using_PcTable
        nset2 = 1;
    else
        nset2 = 3;
    end
end

% Process up to three unique event sets
%   nset = 1 => All events
%   nset = 2 => Screening volume incursion (SVI) events
%   nset = 3 => NonSVI events

for nset=1:nset2
    
    if nset == 1
        eset = true(size(out.ConjDist.Pclud));
        setstr = 'Number of unique events';
        plt_file_root0 = 'All';
    elseif nset == 2
        eset = out.ConjDist.NLTin > 0;
        setstr = 'SVI events';
        plt_file_root0 = 'SVI';
    else
        eset = out.ConjDist.NLTin == 0;
        setstr = 'Non-SVI events';
        plt_file_root0 = 'NonSVI';
    end

    % Number of all events in this set
    Nall = sum(eset);
    setstr = [setstr ': ' num2str(Nall)];
    titstr = setstr;
    
    % Perform processing specifically for the comprehensive set 'All'
    
    if (nset == 1)
        % Last-update option can have RMM rates calculatd
        if PcLUDvsPcMaxOption > 0
            if params.exclude_noncatastrophic
                excatstr = '/catastrophic';
            else
                excatstr = '';
            end
            clear titstrnew;
            titstrnew{1} = titstr;
            if PcLUDvsPcMaxOption == 1
                titstrnew{2} = ['(within RMM commit/consider' excatstr ...
                    ' time limits: ' num2str(sum(out.ConjDist.USElud)) ')'];
            elseif PcLUDvsPcMaxOption == 2
                titstrnew{2} = ['(within RMM commit/consider' excatstr ...
                    ' time limits: ' num2str(sum(out.ConjDist.USEmud)) ')'];
            else
                titstrnew{2} = ' ';
            end
            titstr = titstrnew;
        end
    else
        [p0,dp] = ConjDist_binofit(Nall,Nevent);
        [pmd,plo,phi] = smart_error_range(100*p0,100*dp(1),100*dp(2));
        setstr = [setstr '/' num2str(Nevent)];
        setstr = [setstr ' = ' pmd '%  (95% conf: ' plo '% to ' phi '%)'];
        titstr = [titstr '/' num2str(Nevent) ' = ' pmd ...
            '%  (' plo '% to ' phi '%)'];
        if (nset == 2)
            out.ConjDist.fracs.svi.md = p0;
            out.ConjDist.fracs.svi.lo = dp(1);
            out.ConjDist.fracs.svi.hi = dp(2);
        else
            out.ConjDist.fracs.nosvi.md = p0;
            out.ConjDist.fracs.nosvi.lo = dp(1);
            out.ConjDist.fracs.nosvi.hi = dp(2);
        end
    end        

    logstr = [current_timestring() ' ' setstr];
    log_string(params.logfid,logstr,params.logging,params.displaying);
    
    % Median TCA for events in this set
    Tcset = out.ConjDist.TCAmed(eset);
    
    % Pc valaues for events in this set
    switch params.redyel_event_option
        case 'max'
            Pcset = out.ConjDist.Pcmax(eset);
        case 'min'
            Pcset = out.ConjDist.Pcmin(eset);
        case 'med'
            Pcset = out.ConjDist.Pcmed(eset);
        case 'fud'
            Pcset = out.ConjDist.Pcfud(eset);
        case 'lud'
            Pcset = out.ConjDist.Pclud(eset);
        case 'mud'
            Pcset = out.ConjDist.Pcmud(eset);
        otherwise
            error('Invalid redyel_event_option parameter');
    end
    
    yel = Pcset >= params.Pcyel;
    Nyel = sum(yel);
    [myel,ryel] = ConjDist_binofit(Nyel,Nall);

    red = Pcset >= params.Pcred;
    Nred = sum(red);
    [mred,rred] = ConjDist_binofit(Nred,Nall);
 
    % Make altitude-latitude distribution plots
    
    if (params.make_altlat_plots(nset) ~= 0) && any(eset)

        % Set-specific processing
        if nset == 1
            
            % Initialize plot ranges and name root
            xrng0 = []; yrng0 = []; zrng0 = [];
            plt_file_root = 'All';
            
            % Crate plot title
            nti = 0; titlPlot = cell(1,100);
            nti = nti+1;
            if Npri == 1
                titlPlot{nti} = params.prisetstr;
            else
                titlPlot{nti} = strrep(params.prisetstr,'_',' ');
            end
            
            titlPlot{nti} = ['Orbit: ' titlPlot{nti} ' (' ...
                num2str(round(mdper)) 'km x ' ...
                num2str(round(mdapo)) 'km x ' ...
                num2str(round(mdinc)) 'deg'];

            if isfield(params,'priorb')
                if ~isempty(params.priorb)
                    titlPlot{nti} = [titlPlot{nti} ' in ' params.priorb ')'];
                end
            else
                titlPlot{nti} = [titlPlot{nti} ')'];
            end
            nti=nti+1; titlPlot{nti} = ...
                ['Interval: ' minTCA ' to ' maxTCA ...
                ' (' num2str(delTCA_yr,'%0.2f') ' yr)'];

            nti=nti+1; titlPlot{nti} = ['Commit/consider time limits: ' ...
                num2str(min(params.CommitConsiderDays)) ' to ' ...
                num2str(max(params.CommitConsiderDays)) ' day'];

            if PcLUDvsPcMaxOption > 0
                nti = nti+1;
                if params.exclude_noncatastrophic
                    titlPlot{nti} = ['Likely non-catastrophic events: Excluded (M = ' ...
                        smart_exp_format(params.mission_on_orbit_mass_kg,3) ' kg)'];
                else
                    titlPlot{nti} = 'Likely non-catastrophic events: Included';
                end
            end

            if PcLUDvsPcMaxOption == 1
                nti=nti+1; titlPlot{nti} = ['Unique events: ' ...
                    num2str(sum(out.ConjDist.USElud)) ' of ' num2str(Nall)];
            elseif PcLUDvsPcMaxOption == 2
                nti=nti+1; titlPlot{nti} = ['Unique events: ' ...
                    num2str(sum(out.ConjDist.USEmud)) ' of ' num2str(Nall)];
            else
                nti=nti+1; titlPlot{nti} = ' ';
            end

            if PcLUDvsPcMaxOption > 0
                if params.exclude_noncatastrophic
                    titlPlot{nti} = [titlPlot{nti} ...
                        ' in commit/catastrophic limits'];
                else
                    titlPlot{nti} = [titlPlot{nti} ...
                        ' within commit time limits'];
                end
            end
            altlat_title = titlPlot(1:nti);
            
        elseif nset == 2
            plt_file_root = 'SVI';
        else
            plt_file_root = 'NonSVI';
        end

        % Plot the alt-lat dist in all required plotting modes

        malp = params.make_altlat_plots(nset);
        if malp > 0
            alm1 = 0;
            alm2 = min(3,malp-1);
        else
            alm1 = min(3,abs(malp)-1);
            alm2 = alm1;
        end
        
        for altlat_mode=alm1:alm2
            
            [xrng0,yrng0,zrng0,fighandle2] = ...
                ConjDist_altlat_plot(altlat_mode,PcLUDvsPcMaxOption,eset,altlat_title, ...
                xrng0,yrng0,zrng0,out.ConjDist,params,fighandle2);

            zstr = num2str(altlat_mode);
            
            if (PcLUDvsPcMaxOption > 0) && ...
               (nset == 1) && (altlat_mode == 3)
                pfroot = [output_file_root0 ...
                  '_Dur' num2str(delTCA_yr,'%0.2f') 'yr'];
                plt_file_name = [pfroot '_SpatialDist'];
                plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);
            else
                pfroot = [plt_file_root '_AltLat' zstr '_' ...
                    params.prisetstr '_Ne' num2str(Nall)];
                if params.timestamp
                    pfroot = [pfroot '_' params.timetag];
                end
                plt_file_name = fullfile(params.outputdir,[pfroot params.plot_format]);
            end
            
            saveas(fighandle2,plt_file_name);
            
        end

    end    
    
    % Plot LUD age and Pc distribution histograms
    
    if (nset < 3)
        
        if ishandle(fighandle0)
            set(0, 'CurrentFigure', fighandle0);
        else
            if params.visible_plots
                vis = 'On';
            else
                vis = 'Off';
            end
            fighandle0 = figure_set_up(8.5/11, vis);
        end
        
        % Plot LUD age distribution histogram

        clf;

        AGEdelta = 0.25;
        AGEedges = (0:AGEdelta:max(params.CommitConsiderDays));
        
        if PcLUDvsPcMaxOption == 1
            USEset = out.ConjDist.USElud(eset);
            AGEset = out.ConjDist.AGElud(eset);
        else
            USEset = out.ConjDist.USEmud(eset);
            AGEset = out.ConjDist.AGEmud(eset);
        end
        
        USEable = (USEset == 1);
                
        for plt0=1:-1:0
        
            clf;
            lgnd = cell(4,1);
            hismax = 0;
            
            for icol=0:3

                nlg = icol+1;

                if icol == 0
                    if plt0
                        ndx = USEable;
                    else
                        ndx = USEable & (Pcset >= params.Pcgre);
                    end
                    fcol = params.rcol; ecol = 0.5*fcol;
                    idx = ndx & (Pcset >= params.Pcred);
                    lgnd{nlg} = ['Pc \geq ' smart_exp_format(params.Pcred) ...
                        ' (' num2str(sum(idx)) ' events)'];
                elseif icol == 1
                    ndx = ndx & (Pcset < params.Pcred);
                    fcol = params.ycol; ecol = 0.5*fcol;
                    idx = ndx & (Pcset >= params.Pcyel);
                    lgnd{nlg} = [smart_exp_format(params.Pcyel) ...
                        ' \leq Pc < ' smart_exp_format(params.Pcred) ...
                        ' (' num2str(sum(idx)) ' events)'];
                elseif icol == 2
                    ndx = ndx & (Pcset < params.Pcyel);
                    fcol = params.gcol; ecol = 0.5*fcol;
                    idx = ndx & (Pcset >= params.Pcgre);
                    lgnd{nlg} = [smart_exp_format(params.Pcgre) ...
                        ' \leq Pc < ' smart_exp_format(params.Pcyel) ...
                        ' (' num2str(sum(idx)) ' events)'];
                else
                    ndx = USEable & (Pcset < params.Pcgre);
                    fcol = (2/3)*[1 1 1]; ecol = 0.5*fcol;
                    if plt0
                        lgnd{nlg} = ['Pc < ' smart_exp_format(params.Pcgre) ...
                            ' (' num2str(sum(ndx)) ' events)'];
                    else
                        lgnd{nlg} = ['Pc < ' smart_exp_format(params.Pcgre) ...
                            ' (' num2str(sum(ndx)) ' events, not plotted)'];
                        ndx = [];
                    end
                end

                his = histogram(AGEset(ndx),AGEedges, ...
                    'EdgeColor',ecol,'EdgeAlpha',1, ...
                    'FaceColor',fcol,'FaceAlpha',1);
                hismax = max(hismax,max(his.Values));
                hold on;

            end

            % Plot vertical lines at the remediation maneuver commit time
            % limits

            xrng = [0 max(params.CommitConsiderDays)+AGEdelta];
            yrng = [0 max(hismax+1,round(1.05*hismax))];

            if PcLUDvsPcMaxOption > 0
                lsty = ':';
                lcol = 'k';
                for nn=1:numel(params.CommitConsiderDays)
                    plot([params.CommitConsiderDays(nn) params.CommitConsiderDays(nn)],yrng, ...
                        'LineStyle',lsty,'LineWidth',2, ...
                        'Color',lcol);
                end
            end    

            hold off;

            xlim(xrng); ylim(yrng);
            xlabel('Last Update Age (days)','FontSize',xfsz,'FontWeight',xfwt);
            ylabel('Frequency','FontSize',yfsz,'FontWeight',yfwt);

            clear ttl;
            if Npri == 1
                ttl{1} = params.prisetstr;
            else
                ttl{1} = strrep(params.prisetstr,'_',' ');
            end
            ttl{1} = [ttl{1} ', ' minTCA ' to ' maxTCA];
            ttl{2} = ['HBR = ' num2str(HBRminall)];
            if (HBRminall ~= HBRmaxall)
                ttl{2} = [ttl{2} '-' num2str(HBRmaxall)];
            end
            ttl{2} = [ttl{2} ' m, Commit/Consider Times = ' num2str(min(params.CommitConsiderDays)) ...
                ' to ' num2str(max(params.CommitConsiderDays)) ' days'];
            ttl1 = ttl;
            title(ttl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);

            legend(lgnd,'Location','NorthOutside', ...
                   'FontSize',afsz,'FontWeight',afwt);
            legend('boxoff');

            set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);
            
            if plt0 == 0
                cnjpltstr = '_SIG';
            else
                cnjpltstr = '_ALL';
            end

            if (PcLUDvsPcMaxOption > 0) && (nset == 1)

                plt_file_name = fullfile(params.outputdir, ...
                    [output_file_root '_LastUpdateDist' cnjpltstr ...
                    params.plot_format]);

            else

                plt_file_root = [plt_file_root0 '_AGE_' ...
                    params.prisetstr '_Ne' num2str(Nall) cnjpltstr];
                if params.timestamp
                    plt_file_root = [plt_file_root '_' params.timetag];
                end

                plt_file_name = fullfile(params.outputdir,[plt_file_root params.plot_format]);

            end

            saveas(fighandle0,plt_file_name);
            
        end
        
        % Perform remediation threshold analysis for the specified 
        % prospective mission durations
        
        if (PcLUDvsPcMaxOption > 0) && (nset == 1)
            Tmissionyrs = params.Tmissionyrs;
            Nmissionyrs = numel(Tmissionyrs);
            irp0 = cell(Nmissionyrs,1);
            NPremImage = params.L10PremImage(3);
            L10PremImage = linspace(params.L10PremImage(1), ...
                                    params.L10PremImage(2), ...
                                    NPremImage);
            PremImage = 10.^L10PremImage;
            NTmisImage = params.TmisImage(3);
            TmisImage = linspace(params.TmisImage(1), ...
                                 params.TmisImage(2), ...
                                 NTmisImage);
            L10PcmisTable = params.L10PcmisTable;
            PcmisTable = 10.^L10PcmisTable;
            NPcmisTable = numel(PcmisTable);
            irp = cell(NTmisImage,NPcmisTable);
            TmisTable = params.TmisTable;
            NTmisTable = numel(TmisTable);
        else
            Nmissionyrs = 0;
        end
        
        % Issue an error if Nmissionyrs is not one, which is the expected         
        if (Nmissionyrs ~= 1)
            error('EventRate version of ConjDist requires one single value for Tmissionyrs parameter');
        end
        
        % First-generation remediation threshold analysis and plots
        % (not done by default for EventRate version)
        
        if params.make_first_generation_plots

            % Make a plot of Pcum as a function of Pc

            clf;

            T0yrs = (maxTCAmed-minTCAmed)/365.25;
            Pccut = params.Pc_cutoff_accum_risk;
            cut = (Pcset >= Pccut);

            if PcLUDvsPcMaxOption > 0
                cut = cut & USEset;
            end

            Ncut = sum(cut);
            Tcuse = Tcset(cut);
            Pcuse = Pcset(cut);
            Pcsrt = sort(Pcuse);

            subplot(3,1,1);
            lpc = log10(Pcsrt);
            msiz = 4;
            plot(lpc,(1:Ncut),'x-k','MarkerSize',msiz);
            if PcLUDvsPcMaxOption == 1
                xlabl = 'Log_{10}(Last-update Pc)';
            elseif PcLUDvsPcMaxOption == 2
                xlabl = 'Log_{10}(Max-update Pc)';
            else
                xlabl = ['Log_{10}(Pc,' lower(params.redyel_event_option) ')'];
            end
            xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
            ylabel('Ncum(Tobs)','FontSize',yfsz,'FontWeight',yfwt);
            xrng = [log10(Pccut) 0];
            xlim(xrng);
            yrng = plot_range([1 Ncut],0.05);
            ylim(yrng);

            title(ttl1,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
            set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

            % Calculate bootstrapped cumulative Pc values

            [Boot] = PcCum_PcRem_BootStrap( ...
                T0yrs,Pcuse,Npri,params.Pc_cutoff_accum_risk, ...
                Tmissionyrs,[params.Pcmission params.PcRMMexecution], ...
                params.N_bootstrap_accum_risk,0.05);

            lps = log(1-Pcsrt);
            Prat = 1/Npri;
            lpscum = cumsum(Prat*lps);
            lpccum = log10(1-exp(lpscum));
            LPcCumActualInterp = lpccum(end);

            subplot(3,1,2);
            lpccum0 = lpccum;
            plot(lpc,lpccum0,'x-k','MarkerSize',msiz);
            xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
            ylabel('Log_{10}[Pcum(Tobs)]','FontSize',yfsz,'FontWeight',yfwt);
            xlim(xrng);
            yrng = plot_range(lpccum0);
            ylim(yrng);

            ttl = [];
            [~,Pcumstr,~] = smart_exp_format(10^LPcCumActualInterp,2,[false true]);
            % [~,Pcummdstr,~] = smart_exp_format(10^Boot.LPcCumActualMedian,2,[false true]);
            [~,Pcumlostr,~] = smart_exp_format(10^Boot.LPcCumActualRange(1),2,[false true]);
            [~,Pcumhistr,~] = smart_exp_format(10^Boot.LPcCumActualRange(2),2,[false true]);
            ttl{1} = ['For the Tobs = ' ...
                num2str(T0yrs,'%0.2f') ' yr observation span' ...
                ' Pcum = ' Pcumstr];
            ttl{2} = ['95% resampling range: ' Pcumlostr ...
                ' \leq Pcum \leq ' Pcumhistr];
            title(ttl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
            set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

            Trat = Tmissionyrs/T0yrs;
            lpscum = cumsum(Trat*Prat*lps);
            lpccum = log10(1-exp(lpscum));

            subplot(3,1,3);
            plot(lpc,lpccum,'-k');
            xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
            ylabel('Log_{10}[Pcum(Tm)]','FontSize',yfsz,'FontWeight',yfwt);
            xlim(xrng);

            if PcRemEstimation
                lpccumt = log10(params.Pcmission);
                lpct = interp1(lpccum,lpc,lpccumt,'linear','extrap');
                LPcRemInterp = min(log10(params.Pcmission),lpct);
                disp(['LPcrem: ' num2str(LPcRemInterp) ...
                    ' ' num2str(Boot.LPcRemMedian) ...
                    ' ' num2str(Boot.LPcRemRange(1)) ...
                    ' ' num2str(Boot.LPcRemRange(2))]);
                disp(['Pcrem: ' smart_exp_format(10^LPcRemInterp,3,[false true]) ...
                    ' ' smart_exp_format(10^Boot.LPcRemMedian,3,[false true]) ...
                    ' ' smart_exp_format(10^Boot.LPcRemRange(1),3,[false true]) ...
                    ' ' smart_exp_format(10^Boot.LPcRemRange(2),3,[false true])]);
                yrng = plot_range([min(lpccum) max(lpccum) lpccumt]);
                hold on;
                plot(LPcRemInterp,lpccumt,'s', ...
                    'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
                plot([xrng(1) lpct],[lpccumt lpccumt],'b--');
                plot([lpct lpct],[yrng(1) lpccumt],'b--');
                ylim(yrng);
                ttl = [];
                [~,Pcumstr,~] = smart_exp_format(params.Pcmission,2);
                [~,Premstr,~] = smart_exp_format(10^LPcRemInterp,2,[false true]);
                ttl{1} = ['For a Tm = ' num2str(Tmissionyrs) ' yr mission'];
                ttl{1} = [ttl{1} ' Pcum = ' Pcumstr ...
                    ' with Prem = ' Premstr];
                [~,Premlostr,~] = smart_exp_format(10^Boot.LPcRemRange(1),2,[false true]);
                [~,Premhistr,~] = smart_exp_format(10^Boot.LPcRemRange(2),2,[false true]);
                ttl{2} = ['95% resampling range: ' Premlostr ...
                    ' \leq Prem \leq ' Premhistr];
                ttl2 = [ttl1 ttl];
            else
                lpcremt = log10(params.PcRMMexecution);
                lpct = interp1(lpc,lpccum,lpcremt,'linear','extrap');
                LPcCumInterp = lpct;
                disp(['LPcCum: ' num2str(LPcCumInterp) ...
                    ' ' num2str(Boot.LPcCumMedian) ...
                    ' ' num2str(Boot.LPcCumRange(1)) ...
                    ' ' num2str(Boot.LPcCumRange(2))]);
                disp(['PcCum: ' smart_exp_format(10^LPcCumInterp,3) ...
                    ' ' smart_exp_format(10^Boot.LPcCumMedian,3) ...
                    ' ' smart_exp_format(10^Boot.LPcCumRange(1),3) ...
                    ' ' smart_exp_format(10^Boot.LPcCumRange(2),3)]);
                yrng = plot_range([min(lpccum) max(lpccum) lpct]);
                hold on;
                plot(lpcremt,LPcCumInterp,'s', ...
                    'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
                plot([xrng(1) lpcremt],[lpct lpct],'b--');
                plot([lpcremt lpcremt],[yrng(1) lpct],'b--');
                ylim(yrng);
                ttl = [];
                [~,Pcumstr,~] = smart_exp_format(10^LPcCumInterp,2,[false true]);
                [~,Premstr,~] = smart_exp_format(params.PcRMMexecution,2);
                ttl{1} = ['For a Tm = ' num2str(Tmissionyrs) ' yr mission'];
                ttl{1} = [ttl{1} ' Pcum = ' Pcumstr ...
                    ' with Prem = ' Premstr];
                [~,Pcumlostr,~] = smart_exp_format(10^Boot.LPcCumRange(1),2,[false true]);
                [~,Pcumhistr,~] = smart_exp_format(10^Boot.LPcCumRange(2),2,[false true]);
                ttl{2} = ['95% resampling range: ' Pcumlostr ...
                    ' \leq Pcum \leq ' Pcumhistr];
                ttl2 = [ttl1 ttl];
            end

            title(ttl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
            set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

            clear ttl;

            % Construct file name and save plot
            plt_file_name = [output_file_root '_Pcum1stGen'];
            plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);                
            saveas(fighandle0,plt_file_name);

            % Make plot of PcRem over a range of Tmission values

            if PcRemEstimation && ~isempty(params.Tmissionyrsplot)

                BootOrig = Boot;

                Tma = params.Tmissionyrsplot(1);
                Tmb = params.Tmissionyrsplot(2);
                NTm = params.Tmissionyrsplot(3);
                LTmyrs = linspace(log10(Tma),log10(Tmb),NTm);
                LPcRin = NaN(size(LTmyrs));
                LPcRmd = NaN(size(LTmyrs));
                LPcRlo = NaN(size(LTmyrs));
                LPcRhi = NaN(size(LTmyrs));
                for nTm=1:NTm
                    Tmisyrs = 10^LTmyrs(nTm);
                    [Boot] = PcCum_PcRem_BootStrap( ...
                        T0yrs,Pcuse,Npri,params.Pc_cutoff_accum_risk, ...
                        Tmisyrs,[params.Pcmission params.PcRMMexecution], ...
                        params.N_bootstrap_accum_risk,0.05,BootOrig.ndxBoot);
                    LPcRin(nTm) = Boot.LPcRemInterp;
                    LPcRmd(nTm) = Boot.LPcRemMedian;
                    LPcRlo(nTm) = Boot.LPcRemRange(1);
                    LPcRhi(nTm) = Boot.LPcRemRange(2);
                end

                clf;

                Tmyrs = 10.^LTmyrs;
                loglog(Tmyrs,10.^LPcRlo,':k','LineWidth',1);
                hold on;
                plot(Tmyrs,10.^LPcRhi,':k','LineWidth',1);
                plot(Tmyrs,10.^LPcRin,'-k','LineWidth',1);
                plot(Tmissionyrs,10^LPcRemInterp,'s', ...
                    'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
                plot([Tmissionyrs Tmissionyrs], ...
                     10.^[BootOrig.LPcRemRange(1),BootOrig.LPcRemRange(2)],'b-','LineWidth',2);
                hold off;
                xlabel('Duration (years)','FontSize',xfsz,'FontWeight',xfwt);
                ylabel('RMM Threshold Pc','FontSize',yfsz,'FontWeight',yfwt);

                title(ttl2,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
                set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

                xlim(params.Tmissionyrsplot(1:2));
                yrng = 10.^plot_range([min(LPcRlo) max(LPcRhi)],0.05);
                ylim(yrng);

                plt_file_name = [output_file_root '_Prem1stGen'];
                plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);

                saveas(fighandle0,plt_file_name);

            end

        end % End of make_first_generation_plots

        % Perform risk remediation analysis

        % RMM type can be either translational or rotational
        if strcmpi(params.RemManeuver,'Translational')
            TransMan = true;
        elseif strcmpi(params.RemManeuver,'Rotational')
            TransMan = false;
        else
            error('Invalid RemManuever parameter');
        end

        % Number of years spanned by observed set of conjunctions
        Tobs = (maxTCAmed-minTCAmed)/365.25;

        % Limit usage of observed conjunctions to have Pc above the
        % specified cutoff value
        use = (Pcset >= params.Pc_cutoff_accum_risk);
        if PcLUDvsPcMaxOption > 0
            use = use & USEset;
        end

        % Usable observed conjunctions
        Nobs  = sum(use);
        Pcobs = Pcset(use);

        % Extract mission and processing parameters
        Tmis = Tmissionyrs;
        Nboot = params.N_bootstrap_accum_risk;
        lnPcunc = params.lnPc_bootstrap_uncertainty;
        augevent = params.add_one_unrem_event;
        log10PremSearch = params.log10PremSearch;
        Reduction = params.RemReduction;
        alp95 = (1-0.95)/2;
        alp70 = (1-0.70)/2;
 
        % Plot the unmitigated Pcum as a function of Tmission

        Tma = min(min(params.Tmissionyrs),TmisImage(1));
        Tmb = max(max(params.Tmissionyrs),TmisImage(end));
        NTmisArray = NTmisImage;
        LTmisArray = linspace(log10(Tma),log10(Tmb),NTmisArray);
        TmisArray = 10.^LTmisArray;

        if NTmisTable > 0
            TmisArray = unique([TmisArray TmisTable params.Tmissionyrs]);
        end
        NTmisArray = numel(TmisArray);
        xrng = [TmisArray(1) TmisArray(end)];

        % Calculate unremediated Pcum values
        PcummdNoRem = zeros(NTmisArray,1);
        PcumloNoRem = zeros(NTmisArray,2);
        PcumhiNoRem = zeros(NTmisArray,2);                
        PcRMMval = 1; % No RMMs performed                 
        parfor n=1:NTmisArray % parfor enabled
            Tms = TmisArray(n);
            [pnb] = Prem_Nmvr_Bootstrap(Npri, ...
                Nobs,Tobs,Pcobs,Tms,SecCatGrowth,TransMan,Reduction, ...
                max(1e4,Nboot),[alp70 alp95],lnPcunc,false,PcRMMval,[]);
            PcummdNoRem(n) = pnb.Pcummd;
            PcumloNoRem(n,:) = pnb.Pcumlo';
            PcumhiNoRem(n,:) = pnb.Pcumhi';
        end

        if PcRemEstimation
            % Only make unremediated Pcum plot when estimating Prem
            Nplotpass = 1;
            yrng = 10.^plot_range( ...
                [min(log10(PcumloNoRem(PcumloNoRem > 0))) ...
                 max(log10(PcumhiNoRem))],0.05);
        else
            % Make unremediated and rem. Pcum plots when given Prem
            Nplotpass = 2;
            % Calculate remediated Pcum values
            PcummdWithRem = zeros(NTmisArray,1);
            PcumloWithRem = zeros(NTmisArray,2);
            PcumhiWithRem = zeros(NTmisArray,2);
            PcRMMval = params.PcRMMexecution; % RMMs performed                 
            parfor n=1:NTmisArray % parfor enabled
                Tms = TmisArray(n);
                [pnb] = Prem_Nmvr_Bootstrap(Npri, ...
                    Nobs,Tobs,Pcobs,Tms,SecCatGrowth,TransMan,Reduction, ...
                    max(1e4,Nboot),[alp70 alp95],lnPcunc,false,PcRMMval,[]);
                PcummdWithRem(n) = pnb.Pcummd;
                PcumloWithRem(n,:) = pnb.Pcumlo';
                PcumhiWithRem(n,:) = pnb.Pcumhi';
            end
            yrng = 10.^plot_range( ...
                [min(log10(PcumloNoRem(PcumloNoRem > 0))) ...
                 min(log10(PcumloWithRem(PcumloWithRem > 0))) ...
                 max(log10(PcumhiNoRem)) ...
                 max(log10(PcumhiWithRem))],0.05);
        end

        % Make Pcum Vs Dur plot(s) with one or two passes.
        % (This must be done even if Pcum plots are not being made,
        % to construct the plot title cell array used later for the RMM
        % rate plot.)

        for nplotpass=1:Nplotpass

            if nplotpass == 1
                PcummdPlot = PcummdNoRem;
                PcumloPlot = PcumloNoRem;
                PcumhiPlot = PcumhiNoRem;
            else
                PcummdPlot = PcummdWithRem;
                PcumloPlot = PcumloWithRem;
                PcumhiPlot = PcumhiWithRem;
            end

            if params.make_Pcum_plots
                clf;
                loglog(TmisArray,PcummdPlot,'-k','LineWidth',1);
                hold on;
                % Plot Pcum md,lo,hi curves
                % 95% variation range limits
                plot(TmisArray,PcumloPlot(:,2),':k','LineWidth',1);
                plot(TmisArray,PcumhiPlot(:,2),':k','LineWidth',1);
            end

            % Overplot the mission Pcum estimates

            if ~PcRemEstimation

                % When given RMM red-Pc threshold, estimate the
                % non-remediated Pcum value

                ndx = find(TmisArray == Tmis);
                if ~isempty(ndx)
                    % Mission duration is a member of TmisArray array
                    n = ndx(end);
                    PcummdPL = PcummdPlot(n);
                    PcumloPL = PcumloPlot(n,2);
                    PcumhiPL = PcumhiPlot(n,2);
                else
                    % Mission duration is not a member of TmisArray, so
                    % interpolation must be used
                    ndx = find(TmisArray < Tmis);
                    if isempty(ndx)
                        n1 = 1;
                    else
                        n1 = min(NTmisArray-1,ndx(end));
                    end
                    n2 = n1+1;
                    LTmis = log10(Tmis);
                    LTmisArray = log10(TmisArray((n1:n2))');
                    PcummdPL = 10^interp1(LTmisArray,log10(PcummdPlot(n1:n2))    ,LTmis);
                    PcumloPL = 10^interp1(LTmisArray,log10(PcumloPlot(n1:n2,2))  ,LTmis);
                    PcumhiPL = 10^interp1(LTmisArray,log10(PcumhiPlot(n1:n2,2))  ,LTmis);
                end

                if nplotpass == 1
                    % Mission Pcum estimates with no RMMs
                    PcummdNR = PcummdPL;
                    PcumloNR = PcumloPL;
                    PcumhiNR = PcumhiPL;
                else
                    % Mission Pcum estimates with RMMs
                    PcummdWR = PcummdPL;
                    PcumloWR = PcumloPL;
                    PcumhiWR = PcumhiPL;
                end

                if params.make_Pcum_plots
                    % Plot no-RMM Pcum in blue with vertical error bar
                    plot(Tmis,PcummdPL,'s', ...
                        'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
                    plot([Tmis Tmis],[PcumloPL PcumhiPL],'b-','LineWidth',2);
                end
                
            end

            if params.make_Pcum_plots
                
                hold off;
                % Set axes ranges
                xlim(xrng); ylim(yrng);

                % Overplot grid
                grid on;
                set(gca,'LineWidth',0.35, ...
                    'GridLineStyle','-','GridColor',repmat(0.5,[1 3]), ...
                    'MinorGridLineStyle','-','MinorGridColor',repmat(0.7,[1 3]));

                xlabel('On-Orbit Duration (years)','FontSize',xfsz,'FontWeight',xfwt);
                ylabel('Cumulative Pc','FontSize',yfsz,'FontWeight',yfwt);
                
            end

            % Build titles cell arrays
            
            if nplotpass == 1
                
                % Build title for case with no RMMs performed during
                % mission

                nti = 0; titlPlot = cell(1,100);
                % nti=nti+1; titlPlot{nti} = ' ';
                nti=nti+1; titlPlot{nti} = 'Semi-Empirical Cumulative Collision Probability';
                nti=nti+1; titlPlot{nti} = ' ';

                nti = nti+1;
                if Npri == 1
                    titlPlot{nti} = params.prisetstr;
                else
                    titlPlot{nti} = strrep(params.prisetstr,'_',' ');
                end
                titlPlot{nti} = ['Orbit: ' titlPlot{nti} ' (' ...
                    num2str(round(mdper)) 'km x ' ...
                    num2str(round(mdapo)) 'km x ' ...
                    num2str(round(mdinc)) 'deg'];

                if isfield(params,'priorb')
                    if ~isempty(params.priorb)
                        titlPlot{nti} = [titlPlot{nti} ' in ' params.priorb ')'];
                    end
                else
                    titlPlot{nti} = [titlPlot{nti} ')'];
                end
                nti=nti+1; titlPlot{nti} = ...
                    ['Interval: ' minTCA ' to ' maxTCA ...
                    ' (' num2str(delTCA_yr,'%0.2f') ' yr)'];

                nti=nti+1; titlPlot{nti} = ['Commit/consider time limits: ' ...
                    num2str(min(params.CommitConsiderDays)) ' to ' ...
                    num2str(max(params.CommitConsiderDays)) ' day'];

                if PcLUDvsPcMaxOption > 0
                    nti = nti+1;
                    if params.exclude_noncatastrophic
                        titlPlot{nti} = ['Likely non-catastrophic events: Excluded (M = ' ...
                            smart_exp_format(params.mission_on_orbit_mass_kg,3) ' kg)'];
                    else
                        titlPlot{nti} = 'Likely non-catastrophic events: Included';
                    end
                end

                if PcLUDvsPcMaxOption == 1
                    nti=nti+1; titlPlot{nti} = ['Unique events: ' ...
                        num2str(sum(out.ConjDist.USElud)) ' of ' num2str(Nall)];
                elseif PcLUDvsPcMaxOption == 2
                    nti=nti+1; titlPlot{nti} = ['Unique events: ' ...
                        num2str(sum(out.ConjDist.USEmud)) ' of ' num2str(Nall)];
                else
                    nti=nti+1; titlPlot{nti} = ' ';
                end

                if PcLUDvsPcMaxOption > 0
                    if params.exclude_noncatastrophic
                        titlPlot{nti} = [titlPlot{nti} ...
                            ' in commit/catastrophic limits'];
                    else
                        titlPlot{nti} = [titlPlot{nti} ...
                            ' within commit time limits'];
                    end
                end

                nti=nti+1; titlPlot{nti} = ' ';

                if ~PcRemEstimation
                    nti=nti+1; titlPlot{nti} = [params.mission_name ...
                        ': duration = ' num2str(Tmissionyrs) ...
                        ' years, HBR = ' num2str(params.mission_HBR_meters) ' m' ...
                        ', growth = ' num2str(SecCatGrowth)];
                else
                    nti=nti+1; titlPlot{nti} = ['Mission hard-body radius: ' ...
                        num2str(params.mission_HBR_meters) ' m' ...
                        ', Growth = ' num2str(SecCatGrowth)];
                end
                
                nti=nti+1; titlPlot{nti} = 'No mission risk mitigation maneuvers (RMMs)';

                if ~PcRemEstimation

                    if PcLUDvsPcMaxOption == 1
                        out.CumulativePcNoRMMs = [PcummdNR PcumloNR PcumhiNR];
                    end

                    nti=nti+1; titlPlot{nti} = ' ';

                    [~,Pcumstr,~] = smart_exp_format(PcummdNR,2,[false true]);
                    nti=nti+1; titlPlot{nti} = ['Mission Pc (no RMMs) = ' Pcumstr];
                    [~,Pcumlostr,~] = smart_exp_format(PcumloNR,2,[false true]);
                    [~,Pcumhistr,~] = smart_exp_format(PcumhiNR,2,[false true]);
                    titlPlot{nti} = [titlPlot{nti} ' (95%: ' Pcumlostr ...
                        ' to ' Pcumhistr ')'];
                    txtclr = '{0 0 1}';
                    titlPlot{nti} = ['{\color[rgb]' txtclr titlPlot{nti} '}'];

                end

                nti=nti+1; titlPlot{nti} = ' ';

                % Trim title array
                titlPlot = titlPlot(1:nti);

                PlotStr = 'PcumNoRemVsDur';                        

                if params.make_Pcum_plots
                    logstr = [current_timestring() ...
                        ' Estimated cumulative Pc values with no mission RMMs:'];
                    log_string(params.logfid,logstr,params.logging,params.displaying);
                end

            else
                
                % Output mission Pc with RMMs
                
                if PcLUDvsPcMaxOption == 1
                    out.CumulativePcWithRMMs = [PcummdWR PcumloWR PcumhiWR];
                end
                
                % Build titles for cases of serious-events or with RMMs
                % performed during mission

                % Replace Pcum (no RMMs) with Pcum (with RMMs)
                [~,PremExecutionstr,~] = smart_exp_format(10^L10PremExecution,3);
                nti = 10; titlPlot{nti} = ...
                    ['RMM threshold Pc = ' PremExecutionstr...
                    ' (type = ' lower(params.RemManeuver) ...
                    ', \rho = ' num2str(params.RemReduction) ')'];
                nti=nti+1; titlPlot{nti} = ' ';

                [~,Pcumstr,~] = smart_exp_format(PcummdWR,2,[false true]);
                nti=nti+1; titlPlot{nti} = ['Mission Pc (with RMMs) = ' Pcumstr];
                [~,Pcumlostr,~] = smart_exp_format(PcumloWR,2,[false true]);
                [~,Pcumhistr,~] = smart_exp_format(PcumhiWR,2,[false true]);
                titlPlot{nti} = [titlPlot{nti} ' (95%: ' Pcumlostr ...
                    ' to ' Pcumhistr ')'];
                txtclr = '{0 0 1}';
                titlPlot{nti} = ['{\color[rgb]' txtclr titlPlot{nti} '}'];                
                
                % Change the title with RMM remediation (for later use)
                titlWiRem = titlPlot;
                nti = 10; titlWiRem{nti} = ...
                    ['Serious event threshold Pc = ' PremExecutionstr];

                PlotStr = 'PcumWithRemVsDur';

                if params.make_Pcum_plots
                    logstr = [current_timestring() ...
                        ' Estimated cumulative Pc values with mission RMMs:'];
                    log_string(params.logfid,logstr,params.logging,params.displaying);
                end

            end
            
            if params.make_Pcum_plots

                % Render title
                title(titlPlot,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);

                % Axis fonts
                set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

                % Construct plot file name and save
                plt_file_name = [output_file_root '_' PlotStr];
                plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);
                saveas(gcf,plt_file_name);

                % Report Pcum values
                logstr = ' Duration(yrs)   Pcum(median) Pcum(95%low) Pcum(95%high)';
                log_string(params.logfid,logstr,params.logging,params.displaying);                
                for n=1:NTmisTable
                    ndx = TmisTable(n) == TmisArray;
                    logstr = [                              '    ' ...
                        num2str(TmisTable(n),     '%0.3e')  '    ' ...
                        num2str(PcummdPlot(ndx),  '%0.3e')  '    ' ...
                        num2str(PcumloPlot(ndx,2),'%0.3e')  '    ' ...
                        num2str(PcumhiPlot(ndx,2),'%0.3e')];
                    log_string(params.logfid,logstr,params.logging,params.displaying);
                end
                
            end

        end
        
        % Define the Prem search array (Prem is the RMM Pc threshold)
        L10PremSearch = linspace(log10PremSearch(1), ...
                                 log10PremSearch(2), ...
                                 log10PremSearch(3));
        if ~PcRemEstimation
            % If give RMM Pc, then add it to the search array
            L10PremSearch = unique([L10PremSearch ...
                                    L10PremExecution]);
        end
        % Nsearch = numel(L10PremSearch);
        PremSearch = 10.^L10PremSearch;

        % Plot x-axis range and label
        xrng = [L10PremSearch(1) L10PremSearch(end)];
        xlabl = 'Log_{10}(Threshold Pc)';

        % Calculate the Prem and Nmvr for a bootstrapped array of
        % mission realizations
        [pnb] = Prem_Nmvr_Bootstrap(Npri, ...
            Nobs,Tobs,Pcobs,Tmis,SecCatGrowth,TransMan,Reduction, ...
            Nboot,alp95,lnPcunc,augevent,PremSearch,[]);
        rndboot = pnb.rndboot;
        Pcummd = pnb.Pcummd;
        Pcumlo = pnb.Pcumlo;
        Pcumhi = pnb.Pcumhi;
        Nmvrmd = pnb.Nmvrmd;
        Nmvrlo = pnb.Nmvrlo;
        Nmvrhi = pnb.Nmvrhi;            

        % Make the plot of y = Pcum vs x = Prem

        if PcRemEstimation

            % This branch for PcRem estimation mode (i.e., given
            % mission cumulative Pcum goal, find PcRem and RemRate)

            L10Pcummd = log10(Pcummd);
            L10Pcumlo = log10(Pcumlo);
            L10Pcumhi = log10(Pcumhi);

            irp0 = interp_risk_params(L10Pcmission,L10PremSearch, ...
                                     L10Pcummd,L10Pcumlo,L10Pcumhi, ...
                                     Nmvrmd,Nmvrlo,Nmvrhi);

            L10Premmd = irp0.L10Premmd;
            L10Premlo = irp0.L10Premlo;
            L10Premhi = irp0.L10Premhi;

            IntNmvrmd = irp0.NmvrmdAchieved;
            IntNmvrlo = irp0.NmvrloAchieved;
            IntNmvrhi = irp0.NmvrhiAchieved;

            L10PcummdAchieved = irp0.L10PcummdAchieved;
            L10PcumloAchieved = irp0.L10PcumloAchieved;
            L10PcumhiAchieved = irp0.L10PcumhiAchieved;

            PcumGoalAchieved = ...
                abs(L10PcummdAchieved-L10Pcmission) < 1e-3*L10Pcmission;

            % Make the Pcum plot

            clf;

            % Plot searched solution curves in black
            plot(L10PremSearch,L10Pcummd,'k-','LineWidth',1);
            hold on;                
            plot(L10PremSearch,L10Pcumlo,'k:','LineWidth',1);
            plot(L10PremSearch,L10Pcumhi,'k:','LineWidth',1);

            % Plot the achieved mission-specific solution in blue
            plot(L10Premmd,L10PcummdAchieved,'s', ...
                'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
            plot([L10Premhi L10Premlo],[L10PcummdAchieved L10PcummdAchieved],'b-','LineWidth',2);
            hold off;

            % Axis labels and ranges
            if PcLUDvsPcMaxOption > 0
                ylabl = 'Log_{10}(Pcum)';
            else
                ylabl = ['Log_{10}(Pcum,' lower(params.redyel_event_option) ')'];
            end
            xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
            ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);

            yrng = plot_range([L10Pcumlo L10Pcumhi L10Pcmission],0.05);
            xlim(xrng); ylim(yrng);

            % Build title for PcRemEstimation case
            clear ttl;
            if Npri == 1
                ttl{1} = params.prisetstr;
            else
                ttl{1} = strrep(params.prisetstr,'_',' ');
            end
            ttl{1} = [ttl{1} ', ' minTCA ' to ' maxTCA];

            ttl{2} = ['HBR = ' num2str(HBRminall)];
            if (HBRminall ~= HBRmaxall)
                ttl{2} = [ttl{2} '-' num2str(HBRmaxall)];
            end
            ttl{2} = [ttl{2} ' m, Commit/Consider Times = ' num2str(min(params.CommitConsiderDays)) ...
                ' to ' num2str(max(params.CommitConsiderDays)) ' days'];

            ttl{3} = ['RMM type: ' lower(params.RemManeuver) ...
                ' (\rho = ' num2str(params.RemReduction) ')'];

            [~,Pcmissionstr,~] = smart_exp_format(10^L10Pcmission,2,[false true]);        
            ttl{4} = ['Tmis = ' num2str(Tmissionyrs) 'yr, '...
                      'Pcum goal = ' Pcmissionstr];

            [~,PcummdAchievedstr,~] = smart_exp_format(10^L10PcummdAchieved,2,[false true]);
            [~,PcumloAchievedstr,~] = smart_exp_format(10^L10PcumloAchieved,2,[false true]);
            [~,PcumhiAchievedstr,~] = smart_exp_format(10^L10PcumhiAchieved,2,[false true]);
            ttl{5} = ['Achieved Pcum = ' ...
                PcummdAchievedstr ' (95%: ' PcumloAchievedstr ' to ' ...
                PcumhiAchievedstr ')'];

            [~,Premmdstr,~] = smart_exp_format(10^L10Premmd,2,[false true]);
            [~,Premlostr,~] = smart_exp_format(10^L10Premlo,2,[false true]);
            [~,Premhistr,~] = smart_exp_format(10^L10Premhi,2,[false true]);
            ttl{6} = ['RMM Pc = ' Premmdstr ' (95%: ' Premlostr ...
                ' to ' Premhistr ')'];

            Nremmdstr = num2str(IntNmvrmd/Tmis,'%0.2f');
            Nremlostr = num2str(IntNmvrlo/Tmis,'%0.2f');
            Nremhistr = num2str(IntNmvrhi/Tmis,'%0.2f');
            ttl{7} = ['Estimated RMM Rate = ' Nremmdstr ' / year (95%: ' Nremlostr ...
                ' to ' Nremhistr ')'];

            title(ttl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
            set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

            drawnow;

            % Construct file name and save plot
            plt_file_name = [output_file_root '_PcumVsPrem'];
            plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);
            saveas(fighandle0,plt_file_name);

        end

        % Make the y = RMM-rate vs x = RMM-Pc plot

        clf; plot(NaN,NaN);

        if ~PcRemEstimation
            % Calculated RMM numbers and Pcum values                
            ndx = L10PremSearch == L10PremExecution;
            IntNmvrmd = Nmvrmd(ndx);
            IntNmvrlo = Nmvrlo(ndx);
            IntNmvrhi = Nmvrhi(ndx);
            IntPcummd = Pcummd(ndx);
            IntPcumlo = Pcumlo(ndx);
            IntPcumhi = Pcumhi(ndx);
            [~,PcummdAchievedstr,~] = smart_exp_format(IntPcummd,2,[false true]);
            [~,PcumloAchievedstr,~] = smart_exp_format(IntPcumlo,2,[false true]);
            [~,PcumhiAchievedstr,~] = smart_exp_format(IntPcumhi,2,[false true]);
        end

        % RMM rates
        lorate = Nmvrlo/Tmis;
        mdrate = Nmvrmd/Tmis;
        hirate = Nmvrhi/Tmis;
        IntRatemd = IntNmvrmd/Tmis;

        % Prem x axis is log
        set(gca, 'XScale', 'log');

        % RMM rate y axis linear
        ymx0 = max(hirate);
        ymn0 = min(lorate);
        yrng = plot_range([ymn0 ymx0],0.05); % Linear axis by default
        ylim(yrng);
        IntNmvrloPlot = IntNmvrlo;
        % Change to log y axis for some cases
        yfloor = 0.01;
        if IntNmvrmd >= yfloor
            if ymx0 > 10*ymn0
                ymn1 = floor(log10(max(yfloor,ymn0)));
                ymx1 = log10(max(ymx0,10*ymn1));
                yrng = 10.^[ymn1 ymx1];
                ylim(yrng);
                iidx = lorate <= 0;
                lorate(iidx) = min(yrng(1)/10000,min(lorate(~iidx)));
                iidx = mdrate <= 0;
                mdrate(iidx) = min(yrng(1)/1000,min(mdrate(~iidx)));
                iidx = hirate <= 0;
                hirate(iidx) = min(yrng(1)/100,min(hirate(~iidx)));
                if IntNmvrlo < yrng(1)
                    IntNmvrloPlot = yrng(1);
                end
                set(gca, 'YScale', 'log');
            end
        end
        set(gca, 'YScale', 'log');

        % Plot median and bounding rate curves, which may have been
        % affected by log plotting adjustments
        hold on;
        plot(PremSearch,hirate,'k:','LineWidth',1);
        plot(PremSearch,mdrate,'k-','LineWidth',1);
        plot(PremSearch,lorate,'k:','LineWidth',1);

        if PcRemEstimation

            % Overplot interpolated RMM rate when PcRem is estimated,
            % along with vertical uncertainty bar, both in blue
            plot(10.^L10Premmd,IntRatemd,'s', ...
                'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
            plot(10.^[L10Premmd L10Premmd],[IntNmvrlo IntNmvrhi]/Tmis,'b-','LineWidth',2);

        else
            
            % Estimated mission event rate
            Nremmdval = IntNmvrmd/Tmis;
            Nremloval = IntNmvrlo/Tmis;
            Nremhival = IntNmvrhi/Tmis;
            out.MissionEventRate = [Nremmdval Nremloval Nremhival];

            % Overplot calculated RMM rate when PcRem is provided
            if Nremmdval >= yrng(1)
                NremmdvalPlot = Nremmdval; mrkrPlot = 's';
            else
                NremmdvalPlot = yrng(1); mrkrPlot = 'v';
            end
            plot(10.^L10PremExecution,NremmdvalPlot,mrkrPlot, ...
                'MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1);
            plot(10.^[L10PremExecution L10PremExecution], ...
                [IntNmvrloPlot IntNmvrhi]/Tmis,'b-','LineWidth',2);

            % Generate title for ~PcRemEstimation case
            ttl = titlWiRem;                
            nti = 1; ttl{nti} = 'Semi-Empirical Threshold Pc Event Rates';

            Nremmdstr = num2str(Nremmdval,'%0.2f');
            Nremlostr = num2str(Nremloval,'%0.2f');
            Nremhistr = num2str(Nremhival,'%0.2f');
            nti=12; ttl{nti} = ['Estimated event rate = ' ...
                Nremmdstr ' / year (95%: ' Nremlostr ' to ' Nremhistr ')'];

            txtclr = '{0 0 1}';
            ttl{nti} = ['{\color[rgb]' txtclr ttl{nti} '}'];

            nti=nti+1; ttl{nti} = ' ';

        end

        hold off;

        % Axis limits and labels
        ylabl = 'Event Rate (yr^{-1})';
        if ~(PcLUDvsPcMaxOption > 0)
            ylabl = [ylabl ', ' lower(params.redyel_event_option)];
        end

        xlabel('Threshold Pc','FontSize',xfsz,'FontWeight',xfwt);
        ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
        xlim(10.^xrng);

        % Log grid
        grid on;
        set(gca,'LineWidth',0.35, ...
            'GridLineStyle','-','GridColor',repmat(0.5,[1 3]), ...
            'MinorGridLineStyle','-','MinorGridColor',repmat(0.7,[1 3]));

        % Render title
        title(ttl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);

        % Axis fonts
        set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);
        drawnow

        % Generate file name and save plot
        plt_file_name = [output_file_root '_EventRateVsEventPc'];
        plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);            
        saveas(fighandle0,plt_file_name);
        
        % Write the analysis results
        
        if ~PcRemEstimation

            logstr = ' ';
            log_string(params.logfid,logstr,params.logging,params.displaying);
            logstr = ['--- Summary of Semi-empirical Estimation Analysis ---'];
            log_string(params.logfid,logstr,params.logging,params.displaying);
            
            for iti=2:9
                logstr = titlPlot{iti};
                log_string(params.logfid,logstr,params.logging,params.displaying);
            end

            iti = 10;
            logstr = strrep(titlPlot{iti},'type','RMM type');
            logstr = strrep(logstr,       '\rho','RMM reduction');
            log_string(params.logfid,logstr,params.logging,params.displaying);

            iti = 11;
            logstr = titlPlot{iti};
            log_string(params.logfid,logstr,params.logging,params.displaying);
            
            logstr = ['Mission average above-threshold event rate = ' ...
                Nremmdstr ' / year (95%: ' Nremlostr ' to ' Nremhistr ')'];
            log_string(params.logfid,logstr,params.logging,params.displaying);

            if params.make_Pcum_plots
                
                [~,PcummdAchievedstr,~] = smart_exp_format(IntPcummd,2,[false true]);
                [~,PcumloAchievedstr,~] = smart_exp_format(IntPcumlo,2,[false true]);
                [~,PcumhiAchievedstr,~] = smart_exp_format(IntPcumhi,2,[false true]);
                logstr = ['Mission cumulative Pc (with RMMs) = ' ...
                    PcummdAchievedstr ' (95%: ' PcumloAchievedstr ' to ' ...
                    PcumhiAchievedstr ')'];
                log_string(params.logfid,logstr,params.logging,params.displaying);

                [~,Pcummdstr,~] = smart_exp_format(PcummdNR,2,[false true]);
                [~,Pcumlostr,~] = smart_exp_format(PcumloNR,2,[false true]);
                [~,Pcumhistr,~] = smart_exp_format(PcumhiNR,2,[false true]);
                logstr = ['Mission cumulative Pc (no RMMs) = ' ...
                    Pcummdstr ' (95%: ' Pcumlostr ' to ' Pcumhistr ')'];
                log_string(params.logfid,logstr,params.logging,params.displaying);
                
            end
            
            logstr = ' ';
            log_string(params.logfid,logstr,params.logging,params.displaying);

        end

        % Make Pcum and Nmvr images, and the associated table 

        if params.make_Pcum_images

            % Construct the title for Prem estimation or given Rrem cases

            if PcRemEstimation

                % Construct a title which will also be used for the CSV
                % table
                nti = 0; titltable = cell(1,100);
                nti=nti+1; titltable{nti} = ' ';
                nti=nti+1; titltable{nti} = 'Semi-Empirical Risk Remediation Threshold Analysis';
                nti = nti+1;
                if Npri == 1
                    titltable{nti} = params.prisetstr;
                else
                    titltable{nti} = strrep(params.prisetstr,'_',' ');
                end
                titltable{nti} = ['Orbit: ' titltable{nti} ' (' ...
                    num2str(round(mdper)) 'km x ' ...
                    num2str(round(mdapo)) 'km x ' ...
                    num2str(round(mdinc)) 'deg'];

                if ~isempty(params.priorb)
                    titltable{nti} = [titltable{nti} ' in ' params.priorb ')'];
                else
                    titltable{nti} = [titltable{nti} ')'];
                end
                nti=nti+1; titltable{nti} = ...
                    ['Interval: ' minTCA ' to ' maxTCA ...
                    ' (' num2str(delTCA_yr,'%0.2f') ' yr)'];

                nti=nti+1; titltable{nti} = ['Unique events: ' num2str(Nall)];
                if PcLUDvsPcMaxOption == 1
                    if params.exclude_noncatastrophic
                        titltable{nti} = [titltable{nti} ...
                            ' (' num2str(sum(out.ConjDist.USElud)) ' within commit/catastrophic limits)'];
                    else
                        titltable{nti} = [titltable{nti} ...
                            ' (' num2str(sum(out.ConjDist.USElud)) ' within commit time limits)'];
                    end
                elseif PcLUDvsPcMaxOption == 2
                    if params.exclude_noncatastrophic
                        titltable{nti} = [titltable{nti} ...
                            ' (' num2str(sum(out.ConjDist.USEmud)) ' within commit/catastrophic limits)'];
                    else
                        titltable{nti} = [titltable{nti} ...
                            ' (' num2str(sum(out.ConjDist.USEmud)) ' within commit time limits)'];
                    end
                end

                nti=nti+1; titltable{nti} = ' ';

                nti=nti+1; titltable{nti} = 'Mission & Risk Mitigation Manuever (RMM) Parameters';
                nti=nti+1; titltable{nti} = ['Hard-body protection radius: HBR = ' ...
                    num2str(params.mission_HBR_meters) ' m'];

                nti=nti+1; titltable{nti} = ['RMM type: ' ...
                    lower(params.RemManeuver) ' (reduction factor = ' ...
                    num2str(params.RemReduction) ')'];

                nti=nti+1; titltable{nti} = ['RMM commit time limits: ' ...
                    num2str(min(params.CommitConsiderDays)) ' to ' ...
                    num2str(max(params.CommitConsiderDays)) ' day'];

                nti=nti+1; titltable{nti} = ' ';

                % Trim the table
                titltable = titltable(1:nti);
                titl = titltable;

            else

                % For EventRate mode 1

                titl = titlPlot;                
                nti = 1; titl{nti} = 'Semi-Empirical Cumulative Risk Analysis';
                nti = 11;
                titl(nti) = strrep(ttl(nti+1),'event rate','RMM rate');

                [~,PcummdAchievedstr,~] = smart_exp_format(IntPcummd,2,[false true]);
                [~,PcumloAchievedstr,~] = smart_exp_format(IntPcumlo,2,[false true]);
                [~,PcumhiAchievedstr,~] = smart_exp_format(IntPcumhi,2,[false true]);
                nti=nti+1; titl{nti} = ['Mission Pc (with RMMs) = ' ...
                    PcummdAchievedstr ' (95%: ' PcumloAchievedstr ' to ' ...
                    PcumhiAchievedstr ')'];
                txtclr = '{0 0 1}';
                titl{nti} = ['{\color[rgb]' txtclr titl{nti} '}'];

                [~,Pcummdstr,~] = smart_exp_format(PcummdNR,2,[false true]);
                [~,Pcumlostr,~] = smart_exp_format(PcumloNR,2,[false true]);
                [~,Pcumhistr,~] = smart_exp_format(PcumhiNR,2,[false true]);
                nti=nti+1; titl{nti} = ['Mission Pc (no RMMs) = ' ...
                    Pcummdstr ' (95%: ' Pcumlostr ' to ' Pcumhistr ')'];
                txtclr = '{0 0 1}';
                titl{nti} = ['{\color[rgb]' txtclr titl{nti} '}'];

            end

            % Define plottable quantities

            if PcRemEstimation
                % Estimated RMM Pc threshold values {median & 95% low/high}
                L10PremmdPlt = L10Premmd;
                L10PremloPlt = L10Premmd;
                L10PremhiPlt = L10Premmd;
            else
                % Given RMM Pc threshold value
                L10PremmdPlt = L10PremExecution;
                L10PremloPlt = L10PremExecution;
                L10PremhiPlt = L10PremExecution;
            end

            % Plot L10Fred vs Prem image if required,
            % and only for translational RMMs

            if params.FredVsPremImage && TransMan 

                % Calculate image only for the first mission lifetime

                NL10FredImage = params.L10FredImage(3);
                L10FredImage = linspace(params.L10FredImage(1), ...
                                        params.L10FredImage(2), ...
                                        NL10FredImage);
                PcummdImage0 = NaN(NPremImage,NL10FredImage);
                NmvrmdImage0 = NaN(NPremImage,NL10FredImage);
                parfor n=1:NL10FredImage % parfor enabled
                    Fred = 10^L10FredImage(n);
                    [pnb] = Prem_Nmvr_Bootstrap(Npri, ...
                        Nobs,Tobs,Pcobs,Tmis,SecCatGrowth,TransMan,Fred, ...
                        Nboot,alp95,lnPcunc,augevent,PremImage,rndboot);
                    PcummdImage0(:,n) = pnb.Pcummd';
                    NmvrmdImage0(:,n) = pnb.Nmvrmd' / Tmis;
                end

                ylabl = 'Log_{10}[Reduction Factor \rho]';

                Pcrem_images(L10PremImage,L10PremmdPlt,L10PremloPlt,L10PremhiPlt, ...
                             L10FredImage,log10(Reduction),PcummdImage0, ...
                             [],ylabl,titl);

                % Construct file name and save plot
                plt_file_name = [output_file_root '_RhoVsPremImage'];
                plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);                
                saveas(gcf,plt_file_name);

            end

            % Plot Tmis vs Prem images

            PcummdImage = NaN(NPremImage,NTmisImage);
            NmvrmdImage = NaN(NPremImage,NTmisImage);

            parfor n=1:NTmisImage % parfor enabled

                Tms = TmisImage(n);
                [pnb] = Prem_Nmvr_Bootstrap(Npri, ...
                    Nobs,Tobs,Pcobs,Tms,SecCatGrowth,TransMan,Reduction, ...
                    Nboot,alp95,lnPcunc,augevent,PremImage,[]);

                PcummdImage(:,n) = pnb.Pcummd';

                NmvrmdImage(:,n) = pnb.Nmvrmd' / Tms;

                L10Pcummd = log10(pnb.Pcummd);
                L10Pcumlo = log10(pnb.Pcumlo);
                L10Pcumhi = log10(pnb.Pcumhi);

                for m=1:NPcmisTable
                    irp{n,m} = interp_risk_params( ...
                        L10PcmisTable(m),L10PremImage, ...
                        L10Pcummd,L10Pcumlo,L10Pcumhi, ...
                        pnb.Nmvrmd,pnb.Nmvrlo,pnb.Nmvrhi);
                end

            end

            % Construct version of the images without the nominal
            % solution marked

            if ~PcRemEstimation
                titltable = titl(1:8);
            end

            % Y axis label for the images
            ylabl = 'Duration (years)';

            % Make two-tiered Pcum image (top) RMM-rate image (bottom)
            two_tiered_image = true;
            if two_tiered_image
                % Generate two-tiered image plot
                Pcrem_images(L10PremImage,[],[],[], ...
                             TmisImage,Tmis,PcummdImage,NmvrmdImage, ...
                             ylabl,titltable);
                % Construct file name and save plot
                plt_file_name = [output_file_root '_DurRateVsPremImage'];
                plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);
                saveas(gcf,plt_file_name);
            end

            % Generate single tiered-image of y = Duration vs x = Pcum
            Pcrem_images(L10PremImage,[],[],[], ...
                         TmisImage,Tmis,PcummdImage,[], ...
                         ylabl,titltable);

            % Construct file name and save plot
            plt_file_name = [output_file_root '_DurVsPremImage'];
            plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);
            saveas(gcf,plt_file_name);

            % Revise the title for later use
            if PcRemEstimation
                titltable = [titltable(2:5) titltable(7:end)];
            end

            % Construct a version of the images with the nominal solution
            % marked
            ylabl = 'Duration (years)';
            Pcrem_images(L10PremImage,L10PremmdPlt,L10PremloPlt,L10PremhiPlt, ...
                         TmisImage,Tmis,PcummdImage,[],ylabl,titl);

            % Construct file name and save plot
            plt_file_name = [output_file_root '_DurMisVsPremImage'];
            plt_file_name = fullfile(params.outputdir,[plt_file_name params.plot_format]);            
            saveas(gcf,plt_file_name);

            % Write output tables, if required

            if PcRemEstimation && ...
               (Nmissionyrs > 0) && (NTmisTable > 0) && (NPcmisTable > 0)

                % Write risk remediation tables.
                % Two versions: 1 = nominal estimates
                %               2 = conservative estimates (using 95% ranges)

                for ntabver=1:2

                    if ntabver == 1
                        vstr = '';
                    else
                        vstr = '_95PctUnc';
                    end

                    if nset == 1
                        pfr0 = [ModeStr '_'];
                    else
                        pfr0 = [plt_file_root0 '_'];
                    end

                    if isempty(params.priorb)
                        pos ='';
                    else
                        pos = [strrep(strtrim(params.priorb),' ','_') '_'];
                    end
                    pss = [pos params.prisetstr '_' ...
                        num2str(round(mdper)) 'km' ...
                        num2str(round(mdapo)) 'km' ...
                        num2str(round(mdinc)) 'dg'];

                    tab_file_root = [pfr0 pss ...
                        '_HBR' num2str(params.mission_HBR_meters) 'm' ...
                        '_Commit' num2str(min(params.CommitConsiderDays),'%0.1f') 'dy'];

                    if params.add_one_unrem_event            
                        tab_file_root = [tab_file_root '_Add1ThEv'];
                    end

                    if params.exclude_noncatastrophic
                        tab_file_root = [tab_file_root '_ExcNonCat'];
                    else
                        tab_file_root = [tab_file_root '_IncNonCat'];
                    end

                    tab_file_root = [tab_file_root vstr];

                    tab_file_name = fullfile(params.outputdir,[tab_file_root '.csv']);

                    tabfid = fopen(tab_file_name,'wt');

                    for nti=1:numel(titltable)
                        fprintf(tabfid,'%s\n',titltable{nti});
                    end

                    % Construct Pcum goal list header line

                    outstr = ' , ';
                    for np=1:NPcmisTable
                        outstr = cat(2,outstr, ...
                            [', ' smart_exp_format(PcmisTable(np),3)]);
                        outstr = sprintf([outstr '\t']);
                    end
                    header = outstr;

                    % Tabulate three quantities: Prem, Nmvr and Pcum

                    for nq=1:3

                        if nq == 1
                            outstr = 'Remediation Threshold Pc';
                        elseif nq == 2
                            outstr = 'Remediation Rate (RMM/yr)';
                        else
                            outstr = 'Cumulative Pc Achieved';
                        end
                        fprintf(tabfid,'%s\n',outstr);

                        if nq == 1
                            outstr = ' , , Cumulative Pc Goal';
                            fprintf(tabfid,'%s\n',outstr);                
                            fprintf(tabfid,'%s\n',header);
                        end

                        for nt=1:NTmisTable

                            % Current mission duration for table

                            Tms = TmisTable(nt);

                            % Ensure this table Tmis value is within previously-constructed
                            % image Tmis array

                            ndx = (Tms == TmisImage);
                            Nndx = sum(ndx);
                            if (Nndx == 0)
                                error('All TmisTable values must be in TmisImage array');
                            elseif (Nndx ~= 1)
                                error('Invalid TmisTable or TmisImage value');
                            end

                            % Tabulate data for all mission Pcum goal values

                            if nt == 1
                                outstr = ['Duration (yr) , ' num2str(Tms,'%0.1f')];
                            else
                                outstr = [' , ' num2str(Tms,'%0.1f')];
                            end
                            outstr = sprintf([outstr '\t']);

                            for np=1:NPcmisTable

                                L10Pcm = L10PcmisTable(np);

                                if (L10Pcm == L10Pcmission) && (Tms == Tmis)
                                    % Use interpolated risk parameters for nominal mission
                                    % (which match with the image figure values)
                                    irp1 = irp0;
                                else
                                    % Use interpolated risk parameters saved during image
                                    % construction
                                    irp1 = irp{ndx,np};
                                end

                                tabfmt = '%0.2e';

                                if ntabver == 1
                                    % Nominal remediation recommendations
                                    if nq == 1
                                        qstr = num2str(10^irp1.L10Premmd,tabfmt);
                                    elseif nq == 2
                                        qstr = num2str(irp1.NmvrmdAchieved/Tms,tabfmt);
                                    else
                                        qstr = num2str(10^irp1.L10PcummdAchieved,tabfmt);
                                    end
                                else
                                    % Conservative remediation recommendations
                                    if nq == 1
                                        qstr = num2str(10^irp1.L10Premlo,tabfmt);
                                    elseif nq == 2
                                        qstr = num2str(irp1.NmvrhiAchieved/Tms,tabfmt);
                                    else
                                        qstr = num2str(10^irp1.L10PcumloAchieved,tabfmt);
                                    end
                                end
                                outstr = cat(2,outstr,[', ' qstr]);
                                outstr = sprintf([outstr '\t']);

                            end

                            fprintf(tabfid,'%s\n',outstr);

                        end

                    end

                    % Close CSV table

                    fclose(tabfid);

                end

            end

        end

    end

    % Report SVI rates (typically not done for the EventRate version of
    % ConjDist algorithm)
    
    if params.report_and_plot_svi_rates
    
        % Report SVI rates averaged over entire interval

        out.ConjDist.Nall(nset) = Nall;
        out.ConjDist.Nyel(nset) = Nyel;
        out.ConjDist.Nred(nset) = Nred;

        out.ConjDist.fracs.yel.md(nset) = myel;
        out.ConjDist.fracs.yel.lo(nset) = ryel(1);
        out.ConjDist.fracs.yel.hi(nset) = ryel(2);

        out.ConjDist.fracs.red.md(nset) = mred;
        out.ConjDist.fracs.red.lo(nset) = rred(1);
        out.ConjDist.fracs.red.hi(nset) = rred(2);

        md = Nall;
        sg = chifact*sqrt(Nall);
        out.ConjDist.rates.all.md(nset) = md/out.ConjDist.TCAdelt;
        out.ConjDist.rates.all.lo(nset) = max(0,md-sg)/out.ConjDist.TCAdelt;
        out.ConjDist.rates.all.hi(nset) = max(md+sg)/out.ConjDist.TCAdelt;

        md = Nyel;
        sg = chifact*sqrt(Nyel);
        out.ConjDist.rates.yel.md(nset) = md/out.ConjDist.TCAdelt;
        out.ConjDist.rates.yel.lo(nset) = max(0,md-sg)/out.ConjDist.TCAdelt;
        out.ConjDist.rates.yel.hi(nset) = max(md+sg)/out.ConjDist.TCAdelt;

        md = Nred;
        sg = chifact*sqrt(Nred);
        out.ConjDist.rates.red.md(nset) = md/out.ConjDist.TCAdelt;
        out.ConjDist.rates.red.lo(nset) = max(0,md-sg)/out.ConjDist.TCAdelt;
        out.ConjDist.rates.red.hi(nset) = max(md+sg)/out.ConjDist.TCAdelt;

        [all_md,all_lo,all_hi] = smart_error_range(out.ConjDist.rates.all.md(nset), ...
            out.ConjDist.rates.all.lo(nset),out.ConjDist.rates.all.hi(nset));    
        logstr = ['  Overall average total  rate = ' ...
            all_md ' events/day (95% conf: ' all_lo ' to ' all_hi ')'];
        log_string(params.logfid,logstr,params.logging,params.displaying);

        [yel_md,yel_lo,yel_hi] = smart_error_range(out.ConjDist.rates.yel.md(nset), ...
            out.ConjDist.rates.yel.lo(nset),out.ConjDist.rates.yel.hi(nset));    
        logstr = ['  Overall average yellow rate = ' ...
            yel_md ' events/day (95% conf: ' yel_lo ' to ' yel_hi ')'];
        log_string(params.logfid,logstr,params.logging,params.displaying);

        [red_md,red_lo,red_hi] = smart_error_range(out.ConjDist.rates.red.md(nset), ...
            out.ConjDist.rates.red.lo(nset),out.ConjDist.rates.red.hi(nset));    
        logstr = ['  Overall average red    rate = ' ...
            red_md ' events/day (95% conf: ' red_lo ' to ' red_hi ')'];
        log_string(params.logfid,logstr,params.logging,params.displaying);

        [ymd,ylo,yhi] = smart_error_range(100*myel,100*ryel(1),100*ryel(2));
        logstr = ['  Overall average yellow fraction = ' ...
            ymd '%  (95% conf: ' ylo '% to ' yhi '%)'];
        log_string(params.logfid,logstr,params.logging,params.displaying);

        [rmd,rlo,rhi] = smart_error_range(100*mred,100*rred(1),100*rred(2));
        logstr = ['  Overall average red    fraction = ' ...
            rmd '%  (95% conf: ' rlo '% to ' rhi '%)'];
        log_string(params.logfid,logstr,params.logging,params.displaying);

        % Perform the time binning for this set

        if params.make_time_bins

            % Median TCAs for current subset of events
            TCAset = out.ConjDist.TCAmed(eset);

            % Loop through bins and define rates and fractions
            for nbin=1:out.ConjDist.bins.Nbin

                % Bin time bounds
                t1 = out.ConjDist.bins.Tbinlo(nbin);
                t2 = out.ConjDist.bins.Tbinhi(nbin);
                dt = t2-t1;

                % Find events in current bin, and effective bin time width

                if (nbin < out.ConjDist.bins.Nbin)
                    inbin = (t1 <= TCAset) & (TCAset < t2);
                else
                    inbin = (t1 <= TCAset) & (TCAset <= t2);
                end

                Ninb = sum(inbin);
                md = Ninb;
                sg = chifact*sqrt(Ninb);
                out.ConjDist.rates.binall.md(nset,nbin) = md/dt;
                out.ConjDist.rates.binall.lo(nset,nbin) = max(0,md-sg)/dt;
                out.ConjDist.rates.binall.hi(nset,nbin) = (md+sg)/dt;

                Nyel = sum(inbin & yel);
                md = Nyel;
                sg = chifact*sqrt(Nyel);
                out.ConjDist.rates.binyel.md(nset,nbin) = md/dt;
                out.ConjDist.rates.binyel.lo(nset,nbin) = max(0,md-sg)/dt;
                out.ConjDist.rates.binyel.hi(nset,nbin) = (md+sg)/dt;

                if (Ninb == 0)
                    md = 0; lohi = [0 0];
                else
                    [md,lohi] = ConjDist_binofit(Nyel,Ninb);
                end
                out.ConjDist.fracs.binyel.md(nset,nbin) = md;
                out.ConjDist.fracs.binyel.lo(nset,nbin) = lohi(1);
                out.ConjDist.fracs.binyel.hi(nset,nbin) = lohi(2);

                Nred = sum(inbin & red);
                md = Nred;
                sg = chifact*sqrt(Nred);
                out.ConjDist.rates.binred.md(nset,nbin) = md/dt;
                out.ConjDist.rates.binred.lo(nset,nbin) = max(0,md-sg)/dt;
                out.ConjDist.rates.binred.hi(nset,nbin) = (md+sg)/dt;

                if (Ninb == 0)
                    md = 0; lohi = [0 0];
                else
                    [md,lohi] = ConjDist_binofit(Nred,Ninb);
                end
                out.ConjDist.fracs.binred.md(nset,nbin) = md;
                out.ConjDist.fracs.binred.lo(nset,nbin) = lohi(1);
                out.ConjDist.fracs.binred.hi(nset,nbin) = lohi(2);

            end

            % Estimate bin-to-bin variation median values, and "resample"
            % (i.e., sample using the estimated uncertainties on each point) 
            % to estimate the 95% confidence ranges.
            % Also check if the SVI (nset=2) or nonSVI (nset=3) subsets are
            % identical to the ALL (nset=1) set.  If so, then use those
            % previously resampled results.

            [out.ConjDist.rates.bavgall.md,                 ...
             out.ConjDist.rates.bavgall.lo,                 ...
             out.ConjDist.rates.bavgall.hi] =               ...
                ConjDist_set_resample(nset,           ...
                    out.ConjDist.rates.binall.md,           ...
                    out.ConjDist.rates.binall.lo,           ...
                    out.ConjDist.rates.binall.hi,           ...
                    out.ConjDist.rates.bavgall.md,          ...
                    out.ConjDist.rates.bavgall.lo,          ...
                    out.ConjDist.rates.bavgall.hi,          ...
                    params.bin_samples,0.05,[0 Inf]);

            [out.ConjDist.rates.bavgyel.md,                 ...
             out.ConjDist.rates.bavgyel.lo,                 ...
             out.ConjDist.rates.bavgyel.hi] =               ...
                ConjDist_set_resample(nset,           ...
                    out.ConjDist.rates.binyel.md,           ...
                    out.ConjDist.rates.binyel.lo,           ...
                    out.ConjDist.rates.binyel.hi,           ...
                    out.ConjDist.rates.bavgyel.md,          ...
                    out.ConjDist.rates.bavgyel.lo,          ...
                    out.ConjDist.rates.bavgyel.hi,          ...
                    params.bin_samples,0.05,[0 Inf]);

            [out.ConjDist.rates.bavgred.md,                 ...
             out.ConjDist.rates.bavgred.lo,                 ...
             out.ConjDist.rates.bavgred.hi] =               ...
                ConjDist_set_resample(nset,           ...
                    out.ConjDist.rates.binred.md,           ...
                    out.ConjDist.rates.binred.lo,           ...
                    out.ConjDist.rates.binred.hi,           ...
                    out.ConjDist.rates.bavgred.md,          ...
                    out.ConjDist.rates.bavgred.lo,          ...
                    out.ConjDist.rates.bavgred.hi,          ...
                    params.bin_samples,0.05,[0 Inf]);

            [out.ConjDist.fracs.bavgyel.md,                 ...
             out.ConjDist.fracs.bavgyel.lo,                 ...
             out.ConjDist.fracs.bavgyel.hi] =               ...
                ConjDist_set_resample(nset,           ...
                    out.ConjDist.fracs.binyel.md,           ...
                    out.ConjDist.fracs.binyel.lo,           ...
                    out.ConjDist.fracs.binyel.hi,           ...
                    out.ConjDist.fracs.bavgyel.md,          ...
                    out.ConjDist.fracs.bavgyel.lo,          ...
                    out.ConjDist.fracs.bavgyel.hi,          ...
                    params.bin_samples,0.05,[0 1]);

            [out.ConjDist.fracs.bavgred.md,                 ...
             out.ConjDist.fracs.bavgred.lo,                 ...
             out.ConjDist.fracs.bavgred.hi] =               ...
                ConjDist_set_resample(nset,           ...
                    out.ConjDist.fracs.binred.md,           ...
                    out.ConjDist.fracs.binred.lo,           ...
                    out.ConjDist.fracs.binred.hi,           ...
                    out.ConjDist.fracs.bavgred.md,          ...
                    out.ConjDist.fracs.bavgred.lo,          ...
                    out.ConjDist.fracs.bavgred.hi,          ...
                    params.bin_samples,0.05,[0 1]);

            [all_md,all_lo,all_hi] = smart_error_range(out.ConjDist.rates.bavgall.md(nset), ...
                out.ConjDist.rates.bavgall.lo(nset),out.ConjDist.rates.bavgall.hi(nset));    
            logstr = ['  Bin-to-bin average total  rate = ' ...
                all_md ' events/day (95% conf: ' all_lo ' to ' all_hi ')'];
            log_string(params.logfid,logstr,params.logging,params.displaying);

            [yel_md,yel_lo,yel_hi] = smart_error_range(out.ConjDist.rates.bavgyel.md(nset), ...
                out.ConjDist.rates.bavgyel.lo(nset),out.ConjDist.rates.bavgyel.hi(nset));    
            logstr = ['  Bin-to-bin average yellow rate = ' ...
                yel_md ' events/day (95% conf: ' yel_lo ' to ' yel_hi ')'];
            log_string(params.logfid,logstr,params.logging,params.displaying);

            [red_md,red_lo,red_hi] = smart_error_range(out.ConjDist.rates.bavgred.md(nset), ...
                out.ConjDist.rates.bavgred.lo(nset),out.ConjDist.rates.bavgred.hi(nset));    
            logstr = ['  Bin-to-bin average red    rate = ' ...
                red_md ' events/day (95% conf: ' red_lo ' to ' red_hi ')'];
            log_string(params.logfid,logstr,params.logging,params.displaying);

            [yel_md,yel_lo,yel_hi] = smart_error_range( ...
                100*out.ConjDist.fracs.bavgyel.md(nset), ...
                100*out.ConjDist.fracs.bavgyel.lo(nset), ...
                100*out.ConjDist.fracs.bavgyel.hi(nset));
            logstr = ['  Bin-to-bin average yellow fraction = ' ...
                yel_md '%  (95% conf: ' yel_lo '% to ' yel_hi '%)'];
            log_string(params.logfid,logstr,params.logging,params.displaying);

            [red_md,red_lo,red_hi] = smart_error_range( ...
                100*out.ConjDist.fracs.bavgred.md(nset), ...
                100*out.ConjDist.fracs.bavgred.lo(nset), ...
                100*out.ConjDist.fracs.bavgred.hi(nset));
            logstr = ['  Bin-to-bin average red    fraction = ' ...
                red_md '%  (95% conf: ' red_lo '% to ' red_hi '%)'];
            log_string(params.logfid,logstr,params.logging,params.displaying);

            if params.make_time_plots(nset)

                if nset == 1
                    clear time_title;
                    nti = 1;
                    if Npri == 1
                        time_title{nti} = ['Primary: ' params.prisetstr];
                    else
                        time_title{nti} = ['Primary set: ' params.prisetstr];
                    end
                    time_title{nti} = strrep(time_title{nti},'_',' ');
                    nti=nti+1; time_title{nti} = params.OCMDBroot;
                    time_title{nti} = strrep(time_title{nti},'_','\_');
                    bstr1 = num2str(params.TCAbin_days,'%0.1f');
                    bstr2 = num2str(params.TCAbin_days);
                    if length(bstr1) < length(bstr2)
                        bstr = bstr1;
                    else
                        bstr = bstr2;
                    end
                    nti=nti+1; time_title{nti} = ...
                        ['Span = ' num2str(out.ConjDist.TCAdelt,'%0.1f') ...
                        ' days (Bin = ' bstr ')'];
                    nti=nti+1;       
                    plt_file_root = 'All';
                elseif nset == 2
                    plt_file_root = 'SVI';
                else
                    plt_file_root = 'NonSVI';
                end

                % time_title{nti} = titstr;

                fighandle1 = ConjDist_time_plot(nset,time_title,  ...
                    out.ConjDist.bins,out.ConjDist.rates,out.ConjDist.fracs, ...
                    params,fighandle1);

                plt_file_root = [plt_file_root '_Time_' ...
                    params.prisetstr '_Ne' num2str(Nall)];
                if params.timestamp
                    plt_file_root = [plt_file_root '_' params.timetag];
                end

                plt_file_name = fullfile(params.outputdir,[plt_file_root params.plot_format]);

                % orient(fighandle1,'Landscape');
                saveas(fighandle1,plt_file_name);

            end

        end
        
    end
    
end

%% Close the figures and the log file and return

logstr = [params.start_time ' Function EventRate_ConjDist complete'];
log_string(params.logfid,logstr,params.logging,params.displaying);

if params.logging; fclose(params.logfid); end

return
end
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2022-Apr-11 | Initial version.
% N. Ravago | 2025-Oct-29 | Integrated PcMultiStep
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
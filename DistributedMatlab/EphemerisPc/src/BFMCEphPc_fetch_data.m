function [tca,rca,hca,sca,x1ca,x2ca,Nca,seed,CPUuser,CPUsys,CPUfrac,Hits, ...
    Pday1,Pday2,checkout,nominal,lisTable] = BFMCEphPc_fetch_data(path,mode,Nsamp,ReadOutFiles,ReadSingleFile)
% BFMCEphPc_fetch_data - Get BFMC output data.
%                        (For CARA analysis team internal use)
%
% Syntax: [tca, rca, hca, sca, x1ca, x2ca, Nca, seed, CPUuser, CPUsys, CPUfrac, Hits, ...
%           Pday1, Pday2, checkout, nominal, lisTable] = BFMCEphPc_fetch_data(path, mode, Nsamp, ReadOutFiles, ReadSingleFile);
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
%    path           -   Path to BFMC output directory
%
%    mode           -   Processing mode ('VCM' or 'CDM')
%
%    Nsamp          -   Expected number of MC samples per seed
%
%    ReadOutFiles   -   Flag to read .out files
%
%    ReadSingleFile -   Flag to read a single .lis file 
%
% =========================================================================
%
% Output:
%
%   tca         -   Time of closest approach (s from TCA)
%   rca         -   Relative distance at closest approach
%   hca         -   Hit flags 
%   sca         -   Seed number for each sample
%   x1ca        -   Primary states at CA (m and m/s)
%   x2ca        -   Secondary states at CA (m and m/s)
%   Nca         -   MC sample counts per seed
%   seed        -   Sorted seed values
%   CPUuser     -   User CPU time per seed (s)
%   CPUsys      -   System CPU time per seed (s)
%   CPUfrac     -   CPU fraction per seed
%   Hits        -   Number of hits per seed
%   Pday1       -   Primary propagation time (days)
%   Pday2       -   Secondary propagation time (days)
%   checkout    -   Checkout structure 
%   nominal     -   Nominal conjunction data structure
%   lisTable    -   Concatenated .lis file data table
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Initializations

seed = []; 

tca = []; rca = []; hca =[]; sca = []; x1ca = []; x2ca = []; Nca = []; 
CPUuser = []; CPUsys = []; CPUfrac = [];
Hits = []; Pday1 = []; Pday2 = []; nominal = []; lisTable = [];

checkout.status = -1; checkout.CDMcheckoutSkipped = true;

% VCM mode vs CDM mode

if strcmpi(mode,'VCM')
    VCM_mode = true;
else
    VCM_mode = false;
end

% Augment path
path0 = path;
path = fullfile(path0, [lower(mode) '_results']);    

% Check for a compressed file containing the BFMC outputs

cmpfile = ['bfmc_outputs_' lower(mode) '.tar.gz'];
cmpfile_full = fullfile(path,cmpfile);

cmpfile_used = false;
if exist(cmpfile_full,'file')
    % Extract outputs from compressed file
    disp([' Extracting BFMC outputs from: ' cmpfile]);
    cmpfile_list = untar(cmpfile_full,path);
    cmpfile_used = true;
    disp(['  Number of extracted files: ' num2str(numel(cmpfile_list))]);
end

% Find the BFMC output files

if ReadSingleFile
    fileNumSep = '';
else
    fileNumSep = '_*';
end

template = ['bfmc_' lower(mode) fileNumSep '.lis'];
lisfiles = BFMC_find_data_files(path,template);
Nlisfiles = numel(lisfiles);

if (Nlisfiles == 0)
    lisfiles = BFMC_find_data_files(path0,template);
    Nlisfiles = numel(lisfiles);
    if (Nlisfiles == 0)
        checkout.status = -2;
        warning(['Zero ' template ' files found']);
        return;
    else
        path = path0;
    end
end

% Ensure the first .lis file is readable

f = lisfiles{1};
[~,fff,eee] = fileparts(f);
disp([' Testing validity of data from ' [fff eee]]);
lisData = BFMC_lis_reader(f);

if isempty(lisData)
    checkout.status = -3;
    warning('Invalid or empty .lis file found');
    return;
end

% If the number of samples per .lis file is unknown, then define it
if isempty(Nsamp)
    Nsamp = lisData.MCcounts;
end

% Define the nominal data structure

nominal.PriEpochNum = lisData.PriEpochDN(1);
nominal.PriEpochStr = datestr(nominal.PriEpochNum,'yyyy-mm-dd HH:MM:SS.FFF');
nominal.PriEpochPosVel = [lisData.PriEpochX(1) lisData.PriEpochY(1) lisData.PriEpochZ(1) ...
    lisData.PriEpochXdot(1) lisData.PriEpochYdot(1) lisData.PriEpochZdot(1)]' * 1000;
nominal.PriOrbitPeriod = GetOrbitPeriod(nominal.PriEpochPosVel);
nominal.SecEpochNum = lisData.SecEpochDN(1);
nominal.SecEpochStr = datestr(nominal.SecEpochNum,'yyyy-mm-dd HH:MM:SS.FFF');
nominal.SecEpochPosVel = [lisData.SecEpochX(1) lisData.SecEpochY(1) lisData.SecEpochZ(1) ...
    lisData.SecEpochXdot(1) lisData.SecEpochYdot(1) lisData.SecEpochZdot(1)]' * 1000;
nominal.SecOrbitPeriod = GetOrbitPeriod(nominal.SecEpochPosVel);
nominal.TCAnum = lisData.TCA(1);
nominal.TCAstr = datestr(nominal.TCAnum,'yyyy-mm-dd HH:MM:SS.FFF');
nominal.TCApriposvel = [lisData.PriX(1) lisData.PriY(1) lisData.PriZ(1) ...
    lisData.PriXdot(1) lisData.PriYdot(1) lisData.PriZdot(1)]' * 1000;
nominal.TCAsecposvel = [lisData.SecX(1) lisData.SecY(1) lisData.SecZ(1) ...
    lisData.SecXdot(1) lisData.SecYdot(1) lisData.SecZdot(1)]' * 1000;
nominal.TCAlisposvel = true; % Flag indicating origin from .lis file

% Prepare to read the .out files

if ReadOutFiles

    % Check if any *.out.gz need to be decompressed

    template = ['bfmc_' lower(mode) fileNumSep '.out.gz'];
    outfiles = BFMC_find_data_files(path,template);
    Noutfiles = numel(outfiles);

    if Noutfiles > 0
        for n=1:Noutfiles    
            [ppp,fff,eee] = fileparts(outfiles{n});
            f0 = fullfile(ppp,fff);
            if ~exist(f0,'file') && (n == 1)
                disp([' Decompressing ' fff eee ', etc.']);
            end
            gunzip(outfiles{n});
            delete(outfiles{n});
        end
    end

    % Find the BFMC output files

    template = ['bfmc_' lower(mode) fileNumSep '.out'];
    outfiles = BFMC_find_data_files(path,template);
    Noutfiles = numel(outfiles);

    if (Noutfiles == 0)
        checkout.status = -4;
        warning(['Zero ' template ' files found']);
        return;
    end
    
    if (Nlisfiles ~= Noutfiles)
        error('Unequal numbers of BFMC .lis and .out files');
    end
    
    clear outfiles Noutfiles;
    
end

% Get the CPU time files

template = ['bfmc_' lower(mode) fileNumSep '.tim'];
timfiles = BFMC_find_data_files(path,template);
Ntimfiles = numel(timfiles);

% Ensure that the numbers are the same
if Ntimfiles > 0
    if (Nlisfiles ~= Ntimfiles)
        error('Unequal numbers of BFMC .lis and .tim files');
    end
end

% Number of seeds

Nseed = Nlisfiles;
disp([' Number of unique seeds found = ' num2str(Nseed)]);

% Get the seeds and sort them

seed = NaN(Nseed,1);

for ns=1:Nseed
    [~,fff,~] = fileparts(lisfiles{ns});
    [p,Np] = string_parts(fff,'_');
    seed(ns) = str2double(p(Np));
end

[seed,nsrt] = sort(seed);
lisfiles = lisfiles(nsrt);

% Allocate arrays for samples found per seed and CPU times

Nca     = NaN(Nseed,1);
CPUuser = NaN(Nseed,1);
CPUsys  = NaN(Nseed,1);
CPUfrac = NaN(Nseed,1);
Hits    = NaN(Nseed,1);
Pday1   = NaN(Nseed,1);
Pday2   = NaN(Nseed,1);

% Get the TCA and RCA data for each seed

tca = []; rca = []; sca = []; hca = []; x1ca = []; x2ca = [];

storing_states = []; checkout.status = 0;
bfmc_ver = [];

for ns=1:Nseed
    
    % Seed index
    s  = seed(ns);
    
    % Get data from the .lis files
    f = lisfiles{ns};
    if (ns == 1)
        [~,fff,eee] = fileparts(f);
        disp([' Reading data from ' [fff eee] ', etc.']);
    end
    lisData = BFMC_lis_reader(f);
    if (lisData.MCcounts ~= Nsamp)
        error('Mismatched number of samples');
    end
    
    if ns == 1
        lisTable = lisData;
    else
        lisTable = [lisTable; lisData]; %#ok<AGROW>
    end

    % if VCM_mode
    %     if (lisData.PriPropTime <= 0)
    %         [~,fff,eee] = fileparts(f);
    %         warning(['Nonpositive primary propagation time in ' fff eee]);
    %     end
    %     if (lisData.SecPropTime <= 0)
    %         [~,fff,eee] = fileparts(f);
    %         warning(['Nonpositive secondary propagation time in ' fff eee]);
    %     end
    % end
    
    Nca(ns)   = lisData.MCcounts;
    Hits(ns)  = lisData.MChit;
    
    Pday1(ns) = abs(lisData.PriPropTime);
    Pday2(ns) = abs(lisData.SecPropTime);

    if VCM_mode && (ns == 1)
        
        checkout.status = 0;
  
        checkout.VCMdata = lisData;
        
        [ppp,fff,eee] = fileparts(f);
        fff = strrep(fff,'bfmc_vcm_','bfmc_cdm_');
        if endsWith(ppp,'VCM_BFMC_Data')
            ppp = strrep(ppp,'VCM_BFMC_Data','CDM_BFMC_Data');
        end
        cdmLisFilename = fullfile(ppp,[fff eee]);
        if ~exist(cdmLisFilename,'file')
            disp(['Checkout will be skipped due to missing ' cdmLisFilename]);
            checkout.CDMcheckoutSkipped = true;
        else
            checkout.CDMdata = BFMC_lis_reader(fullfile(ppp,[fff eee]));

            checkout.dTCA     = (checkout.CDMdata.TCA-checkout.VCMdata.TCA)*86400;
            checkout.dPriUsig = checkout.CDMdata.PriUsig-checkout.VCMdata.PriUsig;
            checkout.dPriVsig = checkout.CDMdata.PriVsig-checkout.VCMdata.PriVsig;
            checkout.dPriWsig = checkout.CDMdata.PriWsig-checkout.VCMdata.PriWsig;
            checkout.dSecUsig = checkout.CDMdata.SecUsig-checkout.VCMdata.SecUsig;
            checkout.dSecVsig = checkout.CDMdata.SecVsig-checkout.VCMdata.SecVsig;
            checkout.dSecWsig = checkout.CDMdata.SecWsig-checkout.VCMdata.SecWsig;
            checkout.dRelUsig = checkout.CDMdata.RelUsig-checkout.VCMdata.RelUsig;
            checkout.dRelVsig = checkout.CDMdata.RelVsig-checkout.VCMdata.RelVsig;
            checkout.dRelWsig = checkout.CDMdata.RelWsig-checkout.VCMdata.RelWsig;
            checkout.dRelU    = checkout.CDMdata.RelU   -checkout.VCMdata.RelU   ;
            checkout.dRelV    = checkout.CDMdata.RelV   -checkout.VCMdata.RelV   ;
            checkout.dRelW    = checkout.CDMdata.RelW   -checkout.VCMdata.RelW   ;
            checkout.dRelSep  = checkout.CDMdata.RelSep -checkout.VCMdata.RelSep ;
            checkout.dRelVel  = checkout.CDMdata.RelVel -checkout.VCMdata.RelVel ;
            checkout.CDMcheckoutSkipped = false;
        end
    end
    
    % Get (TCA,RCA) data from the .out files
    
    if ReadOutFiles

        f0 = strrep(f,'.lis','.out');
        if (ns == 1)
            [~,fff,eee] = fileparts(f0);
            disp([' Loading data from ' [fff eee] ', etc.']);
        end
        a = load(f0);
        
        [Nca_ns,Ncol] = size(a);
        
        % Determine the file type that is being read
        
        skipToTimeData = false;
        if isempty(storing_states)
            if (Ncol == 14)
                bfmc_ver = '2.1';
                storing_states = true;
            elseif (Ncol == 2)
                bfmc_ver = '2.1';
                storing_states = false;
            elseif (Ncol == 18)
                bfmc_ver = '3.0';
                storing_states = true;
            elseif (Ncol == 6)
                bfmc_ver = '3.0';
                storing_states = false;
            elseif (Ncol == 0)
                skipToTimeData = true;
            else
                error('Invalid number of columns in .out file');
            end
        end
        
        if ~skipToTimeData
            
            % Record hit flag and number of hits 
            if strcmp(bfmc_ver,'3.0') && ~isempty(a)
                hca = cat(1,hca,a(:,1));
                ind = (a(:,1) == 1);
                Nca_ns = sum(ind);
            end
            
            % Record seed number
            sca = cat(1,sca,repmat(s,[Nca_ns 1]));

            % Record remaining data
            if ~isempty(a)
                if strcmp(bfmc_ver,'2.1')
                    tca = cat(1,tca,a(:,1));
                    rca = cat(1,rca,a(:,2));
                    if storing_states
                        % Convert states to m & m/s
                        x1ca = cat(1,x1ca,a(:,3:8) *1000);
                        x2ca = cat(1,x2ca,a(:,9:14)*1000);
                    end
                else % bfmc_ver == 3.0
                    tca = cat(1,tca,a(:,2));
                    rca = cat(1,rca,a(:,4));
                    if storing_states
                        % Convert states to m & m/s
                        x1ca = cat(1,x1ca,a(:,7:12) *1000);
                        x2ca = cat(1,x2ca,a(:,13:18)*1000);
                    end
                end
            end
            
        end
        
    end
    
    % Get the CPU time data from the .tim files
    
    if Ntimfiles > 0
    
        f0 = strrep(f,'.lis','.tim');
        [lines,Nlines] = get_text_file_lines(f0,1,2);

        uu = NaN; ss = NaN; ff = NaN;

        if Nlines == 2

            % Non-verbose LINUX time command

            [p,Np] = string_parts(lines{1});
            if (Np < 3)
                error('Invalid .tim file contents');
            end
            uu = str2double(strrep(p{1},'user',''));
            ss = str2double(strrep(p{2},'system',''));

            if isnan(uu) || isnan(ss)
                error('Invalid .tim file contents');
            end

        else

            % Verbose LINUX time command

            for nl=1:Nlines
                l = lines{nl};
                k = strfind(l,'User time');
                if ~isempty(k)
                    [p,Np] = string_parts(l);
                    uu = str2double(p{Np});
                else
                    k = strfind(l,'System time');
                    if ~isempty(k)
                        [p,Np] = string_parts(l);
                        ss = str2double(p{Np});
                    else
                        k = strfind(l,'Percent of CPU');
                        if ~isempty(k)
                            [p,Np] = string_parts(l);
                            ff = str2double(strrep(p{Np},'%',''))/100;
                        end                
                    end                
                end
                if ~isnan(uu) && ~isnan(ss) && ~isnan(ff)
                    break;
                end
            end

            if isnan(uu) || isnan(ss) || isnan(ff)
                error('Invalid .tim file contents');
            end

        end

        CPUuser(ns) = uu;
        CPUsys(ns)  = ss;
        CPUfrac(ns) = ff;
        
    end

end

% Ensure that prop times are the same for all of the seeds,
% otherwise issue an error.

if (Nseed > 1)
    
    if any(diff(Pday1) ~= 0)
        error('Primary propagation time mismatch found among .lis files');
    else
        Pday1 = Pday1(1);
    end

    if any(diff(Pday2) ~= 0)
        error('Secondary propagation time mismatch found among .lis files');
    else
        Pday2 = Pday2(1);
    end
    
end

% Delete compressed files

if cmpfile_used
    disp(' Deleting compressed outputs file contents');
    % for nn=1:numel(cmpfile_list); delete(fullfile(path,cmpfile_list{nn})); end
    for nn=1:numel(cmpfile_list); delete(cmpfile_list{nn}); end
end

return
end

function [Period] = GetOrbitPeriod(posvel)
    mu = 3.986004418e5;
    
    rvec = posvel(1:3)' ./ 1000;
    r = norm(rvec);
    vvec = posvel(4:6)' ./ 1000;
    v = norm(vvec);
    
    evec = ((v^2 - mu/r)*rvec - dot(rvec,vvec)*vvec) / mu;
    e = norm(evec);
    
    Eps = v^2/2 - mu/r;
    if e ~= 1.0
        a = -mu/(2*Eps);
        Period = 2*pi*sqrt(a^3/mu);
    else
        error('Cannot calculate orbit period for parabolic orbit!');
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
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
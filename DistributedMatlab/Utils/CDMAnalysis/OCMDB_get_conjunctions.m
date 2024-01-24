function conj = OCMDB_get_conjunctions(params)
%==========================================================================
%
% Search the OCM database and build the primary and secondary states and
% covariances
%
% All output in ECI reference frame and MKS units unless otherwise noted.
%
%==========================================================================

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p, '../CovarianceTransformations')); addpath(s.path);
    s = what(fullfile(p, '../LoggingAndStringReporting')); addpath(s.path);
    pathsAdded = true;
end

% Initializations and defaults

if (nargin == 0); params = []; end

if ~isfield(params,'mat_file'); mat_file = []; else ...
        mat_file = params.mat_file; end
if isempty(mat_file); mat_file = 'OCMDB.mat'; end

if ~isfield(params,'save_UTlimited_mat_file')
    params.save_UTlimited_mat_file = false;
end

if ~isfield(params,'sort_UTlimited_mat_file') || ...
   ~params.save_UTlimited_mat_file
    params.sort_UTlimited_mat_file = false;
end

if ~isfield(params,'UTlimited_created_before_TCA')
    params.UTlimited_created_before_TCA = false;
end

if ~isfield(params,'UTlimited_on_TCA')
    params.UTlimited_on_TCA = true;
end

if ~isfield(params,'UTbegin'); UTbegin = []; else ...
        UTbegin = params.UTbegin; end

if ~isfield(params,'UTend'); UTend = []; else ...
        UTend = params.UTend; end

if ~isfield(params,'primary_set'); primary_set = []; else ...
        primary_set = params.primary_set; end

not_empty_primary_set = ~isempty(primary_set);

if ~isfield(params,'secondary_set'); secondary_set = []; else ...
        secondary_set = params.secondary_set; end

not_empty_secondary_set = ~isempty(secondary_set);

if ~isfield(params,'Pc_20m_cutoff'); Pc_20m_cutoff = []; else ...
        Pc_20m_cutoff = params.Pc_20m_cutoff; end
if isempty(Pc_20m_cutoff); Pc_20m_cutoff = -Inf; end

if ~isfield(params,'DCA_km_cutoff'); DCA_km_cutoff = []; else ...
        DCA_km_cutoff = params.DCA_km_cutoff; end
if isempty(DCA_km_cutoff); DCA_km_cutoff = Inf; end

if ~isfield(params,'Nconj_cutoff'); Nconj_cutoff = []; else ...
        Nconj_cutoff = params.Nconj_cutoff; end
if isempty(Nconj_cutoff); Nconj_cutoff = Inf; end

if ~isfield(params,'full_data_set'); full_data_set = []; else ...
        full_data_set = params.full_data_set; end
if isempty(full_data_set); full_data_set = true; end

if ~isfield(params,'remove_corrupt_data'); remove_corrupt_data = []; else ...
        remove_corrupt_data = params.remove_corrupt_data; end
if isempty(remove_corrupt_data); remove_corrupt_data = false; end

if ~isfield(params,'logfid'); logfid = []; else ...
        logfid = params.logfid; end

if ~isfield(params,'verbose'); verbose = []; else ...
        verbose = params.verbose; end
if isempty(verbose); verbose = 1; end

% Initializations

deg2rad = pi/180;

% Check for a previously saved mat file

[ppp,fff,eee] = fileparts(mat_file);

save_root = fff;

% Check for begin/end date limits

JDfilter = false;
JDbegin = -Inf;
JDend = Inf;

if isempty(UTbegin)
    
    if params.save_UTlimited_mat_file
        error('UTbegin parameter must be specified when saving UT-limited mat file');
    end    
    
else
    
    if length(UTbegin) < 19
        error('Bad date string - needs to be at least 19 characters long');
    end
    
    YEAR    = str2double(UTbegin(1:4));
    MONTH   = str2double(UTbegin(6:7));
    DAY     = str2double(UTbegin(9:10));
    HOUR    = str2double(UTbegin(12:13));
    MINUTE  = str2double(UTbegin(15:16));
    SECOND  = str2double(UTbegin(18:19));
    JDbegin = date2jd(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND);
    
    JDfilter = true;
    
    UTb = UTbegin(3:10);
    UTb = strrep(UTb,'-','');
    
    save_root = [save_root '_Beg' UTb];
    
end

if isempty(UTend)
    
    if params.save_UTlimited_mat_file
        error('UTend parameter must be specified when saving UT-limited mat file');
    end
    
else
    
    if length(UTend) < 19
        error('Bad date string - needs to be at least 19 characters long');
    end
    
    YEAR    = str2double(UTend(1:4));
    MONTH   = str2double(UTend(6:7));
    DAY     = str2double(UTend(9:10));
    HOUR    = str2double(UTend(12:13));
    MINUTE  = str2double(UTend(15:16));
    SECOND  = str2double(UTend(18:19));
    JDend = date2jd(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND);
    
    JDfilter = true;
    
    UTe = UTend(3:10);
    UTe = strrep(UTe,'-','');
    
    save_root = [save_root '_End' UTe];
    
end

if ~isempty(Nconj_cutoff) && (Nconj_cutoff > 0) && ~isinf(Nconj_cutoff)
    save_root = [save_root '_Proc' num2str(Nconj_cutoff)];
else
    save_root = [save_root '_ProcAll'];
end

if ~full_data_set
    save_root = [save_root '_MinData'];
end

save_root = [save_root eee];

save_file = fullfile(ppp,save_root);
if exist(save_file,'file')
    logstr = ['Loading ' save_root];
    log_string(logfid,logstr,[],verbose);    
    load(save_file); % This defines the structure "conj"
    if isfield(conj,'UTbegin') && ...
       isequal(conj.UTbegin,UTbegin) && ...
       isfield(conj,'UTend') && ...
       isequal(conj.UTend,UTend) && ...
       isfield(conj,'Pc_20m_cutoff') && ...
       isequal(conj.Pc_20m_cutoff,Pc_20m_cutoff) && ...       
       isfield(conj,'DCA_km_cutoff') && ...
       isequal(conj.DCA_km_cutoff,DCA_km_cutoff) %#ok<NODEF>
        logstr = ['Loading previously processed and saved OCMDB file: ' save_root];
        log_string(logfid,logstr,[],verbose);       
        return;
    end
end

% Load the OCMDB file file.  This creates a structure called "DB".

logstr = ' ';
log_string(logfid,logstr,[],verbose);

[~,fff,eee] = fileparts(mat_file);    
logstr = ['Loading OCMDB mat file: ' fff eee];
log_string(logfid,logstr,[],verbose);

load(mat_file);

% Ensure proper number of columns

[~,Ncol] = size(DB); %#ok<NODEF>

if Ncol < 218
    error('Expecting at least 218 columns in OCMDB mat file.');
end

% Remove corrupt data if required

ind = DB(:,218) == 1;

if remove_corrupt_data

    DB(ind,:)  = [];
    
    if any(ind)
        logstr = ['  Found and eliminated ' num2str(sum(ind)) ...
            ' corrupt conjunctions of ' num2str(numel(ind)) '.'];
        log_string(logfid,logstr,[],verbose);    
    end
    
else
    
    if any(ind)
        logstr = ['  Found and retained ' num2str(sum(ind)) ...
            ' corrupt conjunctions of ' num2str(numel(ind)) '.'];
        log_string(logfid,logstr,[],verbose);    
    end
    
end

% Filter data out of specified date range

if JDfilter
    
    [NN,~] = size(DB);
    
    ind = false(NN,1);
    
    if params.sort_UTlimited_mat_file
        JDfilt = zeros(NN,1);
    end

    for nn=1:NN
        
        if params.UTlimited_on_TCA
            % OCM TCA date
            JD = date2jd(DB(nn,12), DB(nn,13), DB(nn,14), ...
                         DB(nn,15), DB(nn,16), DB(nn,17));
        else
            % OCM creation date
            JD = date2jd(DB(nn, 3), DB(nn, 4), DB(nn, 5), ...
                         DB(nn, 6), DB(nn, 7), DB(nn, 8));
        end
        
        if params.sort_UTlimited_mat_file
            JDfilt(nn) = JD;
        end
        
        % Mark those out of specified date range, or those not matching
        % specified primary or secondary
        
        if (JD < JDbegin) || (JD > JDend)
            ind(nn) = true;
        elseif not_empty_primary_set && ~ismember(DB(nn,1),primary_set)
            ind(nn) = true;
        elseif not_empty_secondary_set && ~ismember(DB(nn,2),secondary_set)
            ind(nn) = true;
        end
        
    end
    
    % Eliminate those out of date range, and save the raw DB mat file if
    % requested
    
    if any(ind)
    
        % DB(ind,:) = [];
        % JDfilt(ind) = [];
        
        DB = DB(~ind,:);
        JDfilt = JDfilt(~ind,:);
        
        if params.save_UTlimited_mat_file
            
            if params.sort_UTlimited_mat_file
                [~,ind] = sort(JDfilt);
                DB = DB(ind,:);
            end
            
            if params.UTlimited_created_before_TCA
                ind = DB(:,169) > 0;
                DB = DB(ind,:);
            end
        
            filtered_file = fullfile(ppp,[fff '_UTlimited' eee]);
            save(filtered_file,'DB');
            OCMDB_rename_matfile(filtered_file,verbose);
            
        end
        
    end
    
end

% Return now when saving UTlimited files

if params.save_UTlimited_mat_file
    return;
end

% remove data that has Pc (w/ HBR=20m) less than cutoff
ind = DB(:,158) < Pc_20m_cutoff;
DB(ind,:)  = [];

% remove data that has nominal DCA greater than cutoff
ind = DB(:,19) > DCA_km_cutoff;
DB(ind,:)  = [];

% Number of conjunctions

[Nconj,~] = size(DB);

% Display results

logstr = ['  Number of loaded conjunctions = ' num2str(Nconj)];
log_string(logfid,logstr,[],verbose);

% Sort by TCA

logstr = '  Sorting by TCA...';
log_string(logfid,logstr,[],verbose);
[~,ndx] = sort(DB(:,166));
DB = DB(ndx,:);

% logstr = '  Sorting by conjunction creation time...';
% log_string(logfid,logstr,[],verbose);
% [~,ndx] = sort(DB(:,165));
% DB = DB(ndx,:);

% Build the covariance matrix blocks

CovBuild6 = triu(ones(6,6))';
CovBuild6(CovBuild6~=0)  = 1:21;
CovBuild6 = CovBuild6 + triu(CovBuild6',1);

% CovBuild3 = triu(ones(3,3))';
% CovBuild3(CovBuild3~=0)  = 1:6;
% CovBuild3 = CovBuild3 + triu(CovBuild3',1);

% Apply the user-specified maximum number of conjunctions

% Nconj_orig = Nconj;
Nconj = min(Nconj,Nconj_cutoff);

% Allocate the conjunction arrays

logstr = '  Allocating buffers...';
log_string(logfid,logstr,[],verbose);

conj.OCMDB_mat_file = mat_file;

conj.UTC            = cell(1,Nconj);  % UTC at nominal TCA
conj.UTCcreation    = cell(1,Nconj);  % UTC at nominal TCA
conj.id_string      = cell(1,Nconj);  % ID tag of this conjunction

conj.ob1            = NaN(1,Nconj);   % SCN for primary
conj.x1             = NaN(6,Nconj);   % ECI state
conj.C1             = NaN(6,6,Nconj); % ECI covariance

conj.ob2            = NaN(1,Nconj);   % SCN for secondary
conj.x2             = NaN(6,Nconj);   % ECI state
conj.C2             = NaN(6,6,Nconj); % ECI covariance

conj.HBR            = NaN(1,Nconj);   % Combined hard-body radii

if full_data_set
    
    conj.C1_UVW         = NaN(6,6,Nconj); % UVW covariance (e.g., RIC, RSW)
    conj.UTC_ObsLast1   = cell(1,Nconj);  % UTC of last obs
    conj.ObsAgeDays1    = NaN(1,Nconj);   % TCA-T(last obs.)
    conj.apo1           = NaN(1,Nconj);   % Apogee radius
    conj.per1           = NaN(1,Nconj);   % Perigee radius
    conj.inc1           = NaN(1,Nconj);   % Inclination

    conj.C2_UVW         = NaN(6,6,Nconj); % UVW covariance (e.g., RIC, RSW)
    conj.UTC_ObsLast2   = cell(1,Nconj);  % UTC of last obs.
    conj.ObsAgeDays2    = NaN(1,Nconj);   % TCA-T(last obs.)
    conj.apo2           = NaN(1,Nconj);   % Apogee radius
    conj.per2           = NaN(1,Nconj);   % Perigee radius
    conj.inc2           = NaN(1,Nconj);   % Inclination

    conj.Pc_JSPOC       = NaN(1,Nconj);
    conj.Pc_CAS         = NaN(1,Nconj);

    conj.Pc_Foster_20m  = NaN(1,Nconj);
    conj.Pc_Foster      = NaN(1,Nconj);
    conj.Pc_Foster_flag = NaN(1,Nconj);

    conj.Pc_Chan_20m    = NaN(1,Nconj);
    conj.Pc_Chan        = NaN(1,Nconj);
    conj.Pc_Chan_flag   = NaN(1,Nconj);
    
end
   
conj.Corrupt_flag   = NaN(1,Nconj);   % OCMBD corrupt (use/no-use) flag
 
% Loop through conjunctions 

logstr = '  Processing conjunctions...';
log_string(logfid,logstr,[],verbose);

progress_inc  = 10000;
progress_disp = (Nconj > progress_inc);

for nc = 1:Nconj
    
    % Get primary and secondary object numbers
    conj.ob1(nc) = DB(nc,1);
    conj.ob2(nc) = DB(nc,2);
    
    % TCA in UTC, in matlab datenum format
    TCA = DB(nc,166);    
    % TCA in a date string yyyy-mm-dd HH:MM:SS.FFF
    conj.UTC{nc} = datestr(TCA,'yyyy-mm-dd HH:MM:SS.FFF');

    % Tcreation in UTC, in matlab datenum format
    Tcr = DB(nc,165);    
    % TCA in a date string yyyy-mm-dd HH:MM:SS.FFF
    conj.UTCcreation{nc} = datestr(Tcr,'yyyy-mm-dd HH:MM:SS.FFF');
    
    % Combined hard-body radii
    conj.HBR(nc) = DB(nc,157);
    
    if full_data_set
        
        % Primary and secondary last obs.

        conj.UTC_ObsLast1{nc} = datestr(DB(nc,167),'yyyy-mm-dd HH:MM:SS.FFF');
        conj.UTC_ObsLast2{nc} = datestr(DB(nc,168),'yyyy-mm-dd HH:MM:SS.FFF');

        conj.ObsAgeDays1(nc) = DB(nc,170);
        conj.ObsAgeDays2(nc) = DB(nc,171);

        % Primary and secondary orbital params

        conj.apo1(nc) = DB(nc,51) * 1e3;
        conj.per1(nc) = DB(nc,52) * 1e3;
        conj.inc1(nc) = DB(nc,53) * deg2rad;

        conj.apo2(nc) = DB(nc,111) * 1e3;
        conj.per2(nc) = DB(nc,112) * 1e3;
        conj.inc2(nc) = DB(nc,113) * deg2rad;
        
    end
    
    % Primary object ECI state
    X1 = DB(nc,172:177); % [km km/s]
    
    % Primary object covariance 
    TriEls = DB(nc,73:93);
    C1UVW = TriEls(CovBuild6)/1e6; % [km^2 km^2/s km^2/s^2]
    
    % Convert UVW (same as the RIC) to ECI
    C1ECI = RIC2ECI(C1UVW,X1(1:3),X1(4:6));  % [km^2 km^2/s km^2/s^2]
    C1ECI = cov_make_symmetric(C1ECI);
    
    % Secondary object ECI state
    X2 = DB(nc,178:183); % [km km/s]
    
    % Secondary object covariance 
    TriEls = DB(nc,133:153);
    C2UVW = TriEls(CovBuild6)/1e6;  % [km^2 km^2/s km^2/s^2]
    
    % Convert UVW (same as the RIC) to ECI
    C2ECI = RIC2ECI(C2UVW,X2(1:3),X2(4:6));  % [km^2 km^2/s km^2/s^2]
    C2ECI = cov_make_symmetric(C2ECI);
    
    % % Combined covariance in ECI
    % CombCov = C1ECI(1:3,1:3) + C2ECI(1:3,1:3);
    % 
    % % Check the combined covariance against the database covariance
    % TriEls = DB(nc,196:201);
    % CombCovCheck = TriEls(CovBuild3);
    
    % Define the outputs in MKS units
    
    conj.x1(:,nc) = 1e3 * X1;
    conj.x2(:,nc) = 1e3 * X2;
    
    conj.C1(:,:,nc) = 1e6 * C1ECI;    
    conj.C2(:,:,nc) = 1e6 * C2ECI;
    
    if full_data_set
    
        conj.C1_UVW(:,:,nc) = 1e6 * C1UVW;
        conj.C2_UVW(:,:,nc) = 1e6 * C2UVW;

        conj.Pc_JSPOC(nc) = DB(nc,27);

        conj.Pc_CAS(nc) = DB(nc,156);

        conj.Pc_Foster_20m(nc) = DB(nc,158);
        conj.Pc_Foster(nc) = DB(nc,159);
        conj.Pc_Foster_flag(nc) = DB(nc,160);

        conj.Pc_Chan_20m(nc) = DB(nc,161);
        conj.Pc_Chan(nc) = DB(nc,163);
        conj.Pc_Chan_flag(nc) = DB(nc,164);
        
    end
    
    conj.Corrupt_flag(nc) = DB(nc,218);
    
    conj.id_string{nc} = ...
        sprintf('%05i%s%05i%s%04i%02i%02i%s%02i%02i%02i%s%04i%02i%02i%s%02i%02i%02i', ...
        DB(nc,1),'_conj_',DB(nc,2),'_',DB(nc,12:14),'_',DB(nc,15:17),'_',DB(nc,3:5),'_',DB(nc,6:8));

    if progress_disp && (mod(nc,progress_inc) == 0)
        logstr = ['    ' current_timestring() ' processed ' ...
                  num2str(nc) ' of ' num2str(Nconj)];
        log_string(logfid,logstr,[],verbose);        
    end
    
end

if progress_disp
    logstr = ['    ' current_timestring() ' processing complete'];
    log_string(logfid,logstr,[],verbose);
end 

% Save the conj structure for later runs

conj.UTbegin       = UTbegin; 
conj.UTend         = UTend;
conj.Pc_20m_cutoff = Pc_20m_cutoff;
conj.DCA_km_cutoff = DCA_km_cutoff;

logstr = ['Writing file ' save_file];
log_string(logfid,logstr,[],verbose);    

save(save_file,'conj','-v7.3');

logstr = 'Finished writing file';
log_string(logfid,logstr,[],verbose);    

return;
end

function Csym = cov_make_symmetric(C)

% Make a covariance matrix symmetric.

szC = size(C);

if (numel(szC) ~= 2)
    error('Array needs to be a 2D matrix to make symmetric.');
end

if (szC(1) ~= szC(2))
    error('Matrix needs to be square to make symmetric.');
end

if ~issymmetric(C)
    Csym = (C+C')/2;
else
    Csym = C;
end

return;
end
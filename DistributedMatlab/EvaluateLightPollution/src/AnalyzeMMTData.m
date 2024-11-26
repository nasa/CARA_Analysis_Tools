function out = AnalyzeMMTData(datapath,Rnormkm,ExtinctionCoef,UTbegin,UTend,params)
% Analyze MMT Data - Generate CDFs from MMT photometric data.
% Syntax: out = AnalyzeMMTData(datapath,Rnormkm,ExtinctionCoef);
%         out = AnalyzeMMTData(datapath,Rnormkm,ExtinctionCoef,UTbegin);
%         out = AnalyzeMMTData(datapath,Rnormkm,ExtinctionCoef,UTbegin,UTend);
%         out = AnalyzeMMTData(datapath,Rnormkm,ExtinctionCoef,UTbegin,UTend,params);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Fit a three-component model to MMT photometric data, and 
% generate associated CDFs.
%
% =========================================================================
%
% Input:
%
%    datapath       - [String] Path to directory containing MMT data files.
%                   
%    Rnormkm        - [Scalar] Range normalization distance (km).
%
%    ExtinctionCoef - [Scalar] Atmospheric extinction coefficient.
%
%    UTbegin        - [String] (Optional) Start of observational period in 
%                     Universal Time. Defaults to all available data if not
%                     provided.
%
%    UTend          - [String] (Optional) End of observational period in 
%                     Universal Time. Defaults to all available data if not
%                     provided.
%
%    params         - [Structure] Auxiliary parameter structure with
%                     optional settings fields.
%
% =========================================================================
%
% Output:
%
%   out             - [Structure] Output structure containing analysis
%                      results.
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

% Initializations
Nargin = nargin;
if Nargin < 6; params  = []; end
if Nargin < 5; UTend   = ''; end
if Nargin < 4; UTbegin = ''; end

% Default parameters

% Check if previous data has been input
if ~isfield(params,'previous_output'); params.previous_output = []; end
previous_output_exists = ~isempty(params.previous_output);

% Flux scaling factor, usually to applt to previously analyzed data
if ~isfield(params,'Fscale') || isempty(params.Fscale)
    Fscale = 1;
else
    Fscale = params.Fscale;
end

% Statistical analysis flag to estimate quantiles and CDF parameters
if ~isfield(params,'StatAnalysis') || isempty(params.StatAnalysis)
    StatAnalysis = true;
else
    StatAnalysis = params.StatAnalysis;
end

% Number of aA model components; zero for no model
if ~isfield(params,'Ncomp') || isempty(params.Ncomp)
    Ncomp = 3;
else
    Ncomp = params.Ncomp;
    if ~(Ncomp == 0 || Ncomp == 1 || Ncomp == 2 || Ncomp == 3)
        error('Invalid number of model aA fitting components');
    end
end

% Output path
if ~isfield(params,'outpath') || isempty(params.outpath)
    outpath = datapath;
else
    outpath = params.outpath;
end

% Mat file saving flag
if ~isfield(params,'UseMMTSaveFiles') || isempty(params.UseMMTSaveFiles)
    UseMMTSaveFiles = ~previous_output_exists;
else
    UseMMTSaveFiles = params.UseMMTSaveFiles;
end

% Plotting flag
if ~isfield(params,'Plotting') || isempty(params.Plotting)
    Plotting = ~previous_output_exists;
else
    Plotting = params.Plotting;
end
if Plotting
    StatAnalysis = true;
end

% Verbose flag
if ~isfield(params,'verbose') || isempty(params.verbose)
    verbose = ~previous_output_exists;
else
    verbose = params.verbose;
end

% Model components (up to three)
if (Ncomp > 0)
    comp = {'SpecSphere' 'DiffSphere' 'DiffPABNormalFacet'};
    comp = comp(1:Ncomp);
end

% CDF reporting levels
CDFrep = [0.50 0.90 0.95];
ndxstar = 2; % Index for reference zentith magniture 

% Normalization range squared
Rnormm2 = (Rnormkm*1e3)^2;

% Clear filter solar magnitude
Msun1AU = -26.91;
AUkm = 149597870.700; % 1.49598e8;

% JD minus Matlab datenum
JDminusDN = 1721058.5;

% Number of digits for fractions
Ndig = 4;

% Initialized data structure
dat = [];

% Generate output file name root
[~,datpth,epth] = fileparts(datapath);
if ~isempty(epth)
    datpth = [datpth '.' epth];
end
outroot = [datpth ...
           '_' num2str(Rnormkm) 'km' ...
           '_Sc' smart_exp_format(Fscale,Ndig)];

% Check if a previous run has been provided
if previous_output_exists
    
    % Extract MMT output from the previous run
    out = params.previous_output;
    
else
    
    % Search for a previously saved .mat file
    need_to_read_data = true;
    if UseMMTSaveFiles
        dl = dir(fullfile(outpath,[outroot '*_AnalyzeMMTData.mat']));
        Nfile = numel(dl);
        if Nfile == 1
            out = []; load(fullfile(outpath,dl.name));
            if isequal(out.datapath , datapath ) && ...
               isequal(out.Rnormkm  , Rnormkm  ) && ...
               isequal(out.Fscale   , Fscale   ) && ...
               isequal(out.UTbegin  , UTbegin  ) && ...
               isequal(out.UTend    , UTend    ) && ...
               isequal(out.Ncomp    , Ncomp    )
                need_to_read_data = false;
            end
        elseif Nfile ~= 0
            error('More than one *_AnalyzeMMTData.mat file found');
        end
    end

    % Use the mat file if appropriate
    if verbose
        if ~need_to_read_data
            disp(['Using save file ' dl.name]);
        else
            disp(['Processing MMT data in ' datapath]);
        end
    end
    
    if need_to_read_data
        
        % Get all text files in data directory
        dl = dir(fullfile(fullfile(datapath),'satellite*.txt'));
        Nfile = numel(dl);
        % Issue an error if no data files are found
        if Nfile == 0
            error('No photometric data files found (they may need to be unzipped)');
        elseif verbose
            disp(['Reading ' num2str(Nfile) ' data files']);
        end
        
        % Date numbers
        if isempty(UTbegin)
            DNbegin = -Inf;
        else
            DNbegin = datenum(UTbegin,'yyyy-mm-dd HH:MM:SS.FFF');
        end
        if isempty(UTend)
            DNend = Inf;
        else
            DNend   = datenum(UTend,'yyyy-mm-dd HH:MM:SS.FFF');
        end

        % Get all of the data, eliminating out-of-range dates and penumbral/umbral
        % data points
        dat = [];
        for n=1:Nfile
            % Read the MMT data table
            [dt,~,~] = ReadMMTData(fullfile(datapath,dl(n).name));
            % MMT satellite number
            [~,fff,~] = fileparts(dl(n).name);
            [prt,~] = string_parts(fff,'_');
            MMTSat = str2double(prt{2});
            % Only keep desired data
            ndx = (dt.Penumbra == 0)      & ...
                  (dt.DateNum >= DNbegin) & ...
                  (dt.DateNum <= DNend);
            dt = dt(ndx,:);
            Ndt = size(dt,1);
            Ntr = numel(unique(dt.Track));
            if verbose
                disp([' MMTSat = ' num2str(MMTSat) ...
                    ', Ntrk = ' num2str(Ntr) ', Ndat = ' num2str(Ndt)]);
            end
            if Ndt == 0
                break;
            end
            % Add the MMT satellite number to the table   
            dt.MMTSat = repmat(MMTSat,[Ndt 1]);
            % Concatenate the data
            if isempty(dat)
                dat = dt;
            else
                dat = vertcat(dat, dt); %#ok<AGROW>
            end
        end

        % Add solar distances to the data
        Ndat = size(dat,1);
        dat.SunDistance = NaN(Ndat,1);
        for n=1:Ndat
            JD = dat.DateNum(n)+JDminusDN;
            rsun = SunPos(JD);
            dat.SunDistance(n) = norm(rsun);
        end
        
        % Set data Rnormkm and Fscale to placeholder values to indicate
        % that further analysis will be required for this data set
        out.Rnormkm = -1;
        out.Fscale  = -1;
        
    end

end

% Check if any further analysis is required
if ~params.Plotting && (Rnormkm == out.Rnormkm) && (Fscale == out.Fscale)
    % Return the current out structure
    return;
else
    % Extract data for further analysis, if required
    if isempty(dat)
        dat = out.dat;
    end
end

% Calculate the range-normalized mag and OCS, including the flux scale
% factor if required
if verbose
    disp([' Normalizing to a range of ' num2str(Rnormkm) ' km']);
    if Fscale ~= 1
        disp([' Adjusting by a flux scaling factor of ' num2str(Fscale)]);
    end
end
dat.MagNorm = dat.Mag ...
    - 5.0 * log10(dat.Distance/Rnormkm) ...
    - 5.0 * log10(dat.SunDistance/AUkm) ...
    - 2.5 * log10(Fscale);

% Total number of data points, tracks and satellites
Ndat = size(dat,1);
Ntrk = numel(unique(dat.Track));
Nsat = numel(unique(dat.MMTSat));

if verbose
    disp('---- TOTAL DATA SET ----');
    disp(['Nsat = ' num2str(Nsat) ', Ntrk = ' num2str(Ntrk) ', Ndat = ' num2str(Ndat)]);
end

% Define output structure
out.datapath = datapath;
out.Rnormkm = Rnormkm;
out.Fscale = Fscale;
out.UTbegin = UTbegin;
out.UTend = UTend;
out.Ncomp = Ncomp;
out.Nsat = Nsat;
out.Ntrk = Ntrk;
out.Ndat = Ndat;
out.dat = dat;
out.CDFfrac = CDFrep;

% Calculate the magnitude CDF and quantiles
if StatAnalysis
    [F,X] = ecdf(dat.MagNorm);
    F = flipud(1-F);
    NF = numel(F);
    X = flipud(X);
    mrep = NaN(size(CDFrep));
    Nrep = numel(CDFrep);
    for nrep=1:Nrep
        ndx = find(F < CDFrep(nrep));
        ndx = ndx(end);
        if (ndx == NF)
            idx = [ndx-1 ndx];
        else
            idx = [ndx ndx+1];
        end
        mrep(nrep) = interp1(F(idx),X(idx),CDFrep(nrep),'linear','extrap');
    end
    out.CDFMagNorm = mrep;
    out.MagNormStar = mrep(ndxstar);
    % Calculate the +1sigma and -1sigma equivalent deviations to match
    % central 90%
    mmd = mrep(1);
    mlo = mrep(3);
    ndx = find(F < 0.05);
    ndx = ndx(end);
    if (ndx == NF)
        idx = [ndx-1 ndx];
    else
        idx = [ndx ndx+1];
    end
    mhi = interp1(F(idx),X(idx),0.05,'linear','extrap');
    siglo = (mlo-mmd)/sqrt(2)/erfinv(-1+2*0.05);
    sighi = (mhi-mmd)/sqrt(2)/erfinv(-1+2*0.95);
    out.StatsMagNorm = [mmd siglo sighi];
end

% Calculate SATCON-1 cutoff
mcut = SatCon1Limit(Rnormkm);
out.MagnitudeSatCon1Limit = mcut;

% Calculate fraction of zenith mags brighter than cutoff
fcut = sum(dat.MagNorm < mcut)/numel(dat.MagNorm);
out.FracBrighterThanSatCon1Limit = fcut;

% Calculate fraction of off-zenith mags bright than cutoff
Nzag = numel(params.ZenAngGrid);
Rekm = 6378.137;
akm = Rekm+Rnormkm; akm2 = akm^2;
out.ZenAngGrid = params.ZenAngGrid;
out.FracBrightGrid0 = NaN(size(params.ZenAngGrid));
out.FracBrightGrid  = NaN(size(params.ZenAngGrid));
AirMass = AirMassRozenberg1966(params.ZenAngGrid);
for nzag=1:Nzag
    czag = cos(params.ZenAngGrid(nzag)); szag = sin(params.ZenAngGrid(nzag));
    rho = sqrt(akm2-(Rekm*szag)^2) - Rekm*czag;
    % AirMass = 1/czag;
    dM = 5*log10(rho/Rnormkm);
    out.FracBrightGrid0(nzag) = sum(dat.MagNorm + dM < mcut);    
    dM = dM + ExtinctionCoef*AirMass(nzag);
    out.FracBrightGrid(nzag) = sum(dat.MagNorm + dM < mcut);
end
out.FracBrightGrid0 = out.FracBrightGrid0/numel(dat.MagNorm);
out.FracBrightGrid  = out.FracBrightGrid /numel(dat.MagNorm);

% Calculate the multicomponent aA model

if Ncomp > 0

    if verbose
        disp(['Fitting the data using ' num2str(Ncomp) ' component aA model']);
    end
    
    % OCS values
    dat.OCSm2 = Rnormm2 * 10.^(-0.4*(dat.MagNorm-Msun1AU));
    
    % Generate the design matrix aA products
    Dmat = NaN(Ndat,Ncomp);
    for n=1:Ncomp
        Dmat(:,n) = ThreeCompComponent(comp{n},dat);
    end

    % Estimate nonnegative aA products
    aAvec = lsqnonneg(Dmat,dat.OCSm2);

    % Analyze to find lambda = log(aA) values for Least Absolute deviation
    logy = log(dat.OCSm2);
    errfcn = @(lambda)sum(abs(logy-log(Dmat*exp(lambda)))); % L1 error-function 

    fmsopts = optimset('fminsearch');
    fmsopts.MaxIter     = 2000*Ncomp;
    fmsopts.MaxFunEvals = 2*fmsopts.MaxIter;

    aAscale = 1.5*sum(aAvec); aAmin = 1e-3*aAscale;
    Nsearch = 30;
    lsearch = NaN(Ncomp,Nsearch);
    esearch = NaN(1,Nsearch);

    for ns=1:Nsearch
        aA0 = aAscale*rand(Ncomp,1);
        aA0(aA0 < aAmin) = aAmin;
        lam0 = log(aA0);
        [lsearch(:,ns),esearch(ns)] = fminsearch(errfcn,lam0,fmsopts); % L1-optimization
    end
    [~,ns] = min(esearch);
    aAvec = exp(lsearch(:,ns));

    sumaAvec = sum(aAvec);
    ndx = aAvec < 1e-6*sumaAvec;
    if any(ndx)
        aAvec(ndx) = 0;
        ndx = ~ndx;
        errfcn = @(lambda)sum(abs(logy-log(Dmat(:,ndx)*exp(lambda))));
        lam1 = fminsearch(errfcn,lsearch(ndx,ns),fmsopts); % L1-optimization
        aAvec(ndx) = exp(lam1);
    end
    OCSmod = Dmat * aAvec;
    
    out.aA = aAvec;
    
    % Analyze the distibution of OCS value above the model as a function of aA
    % scaling factor

    f1 = 0;
    f2 = 2.5;
    Nf = 1000;
    CDF = NaN(1,Nf);
    maxCDFrep = max(CDFrep);
    iterating = true;
    while iterating
        f = linspace(f1,f2,Nf);
        for n=1:Nf
            ndx = dat.OCSm2 <= f(n)*OCSmod;
            CDF(n) = sum(ndx)/Ndat;
        end
        if CDF(end) < maxCDFrep
            f2 = 2*f2;
        else
            iterating = false;
        end
    end

    frep = NaN(size(CDFrep));
    for nrep=1:Nrep
        ndx = find(CDF < CDFrep(nrep));
        ndx = ndx(end);
        if (ndx == Nf)
            idx = [ndx-1 ndx];
        else
            idx = [ndx ndx+1];
        end
        frep(nrep) = interp1(CDF(idx),f(idx),CDFrep(nrep),'linear','extrap');
    end
    out.CDFaAScaleFact = frep;
    
end

% Generate output file name
outfile = [outroot                                  ...
           '_' datestr(min(dat.DateNum),'yyyymmdd') ...
           '_' datestr(max(dat.DateNum),'yyyymmdd') ...
           '_Ns' num2str(Nsat) '_Nt' num2str(Ntrk)  ...
           '_Nd' num2str(Ndat) '_Nc' num2str(Ncomp)];

% Set up for plotting

if Plotting
    
    figure(1); clf;

    xfsz = 11;
    yfsz = xfsz;

    xfwt = 'bold';
    yfwt = xfwt;

    tfsz = xfsz;
    tfwt = 'bold';
    tfan = 'italic';

    afsz = xfsz;
    afwt = 'bold';
    alwd = 0.5;

    msiz = 2;
    
    % Generate title
    titl = cell(1,1);
    nt = 1;  titl{nt} = [strrep(datpth,'_','\_') ...
        ' (Nsat = ' num2str(Nsat) ', Ntrk = ' num2str(Ntrk) ...
        ', Ndat = ' num2str(Ndat) ')'];
    UTearliest = datestr(min(dat.DateNum),'yyyy-mm-dd HH:MM:SS');
    UTlatest   = datestr(max(dat.DateNum),'yyyy-mm-dd HH:MM:SS');
    nt=nt+1; titl{nt} = ['Earliest obs. date = ' UTearliest];
    nt=nt+1; titl{nt} = ['Latest   obs. date = ' UTlatest];
    if Fscale ~= 1
        nt=nt+1; titl{nt} = ['Flux scale factor = ' ...
            smart_exp_format(Fscale,Ndig)];
    end
    nt0 = numel(titl);
    
    % Set x-axis range
    xrng = plot_range(dat.Phase,0.05);
    xrng(xrng < 0) = 0;
    xrng(xrng > 180) = 180;

    % Plot the multicomponent OCS model, if required    
    if Ncomp > 0
    
        clf;
        mrkr = 'o';
        lwid = 0.5;
        semilogy(dat.Phase,dat.OCSm2, ...
                'LineStyle','none', ...
                'LineWidth',lwid, ...
                'Marker',mrkr, ...
                'MarkerSize',msiz, ...
                'MarkerFaceColor','k', ...
                'MarkerEdgeColor','k');
        hold on;
        [~,srt] = sort(dat.Phase);
        plot(dat.Phase(srt),OCSmod(srt),'-r','LineWidth',1);
        hold off;
        xlim(xrng);
        xlabl = 'Phase Angle (deg)';
        ylabl = 'OCS (m^2)';
        xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
        ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);

        nt=nt+1; titl{nt} = ['Least absolute magnitude deviation ' ...
            num2str(Ncomp) '-component solution'];
        for n=1:Ncomp
            titl{nt+n} = [comp{n} ': aA = ' ...
                smart_exp_format(aAvec(n),3) ' m^2'];
        end
        nt=numel(titl)+1; titl{nt} = ' ';        
        title(titl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
        
        set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

        saveas(gcf,[fullfile(outpath,outfile) '_OCSPhase.png']);

        nt=numel(titl); titl{nt} = ['Scaling factors: ' ...
            smart_exp_format(frep(1),3) ' (' smart_exp_format(100*CDFrep(1),3) '%) ' ...
            smart_exp_format(frep(2),3) ' (' smart_exp_format(100*CDFrep(2),3) '%) ' ...
            smart_exp_format(frep(3),3) ' (' smart_exp_format(100*CDFrep(3),3) '%)'];
        nt=nt+1; titl{nt} = ' ';

        clf;
        plot(f,CDF,'-k','LineWidth',1);
        grid on;
        xlabl = 'aA Scaling Factor';
        ylabl = 'CDF';
        xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
        ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
        title(titl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
        set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

        saveas(gcf,[fullfile(outpath,outfile) '_aACDF.png']);
        
    end

    % Plot a CDF of range-normalized magnitudes
    % Put quantiles in title
    titl = titl(1:nt0); nt = numel(titl);
    nt=nt+1; titl{nt} = ['Magnitude quantiles: ' ...
        smart_exp_format(mrep(1),3) ' (' num2str(100*CDFrep(1)) '%) ' ...
        smart_exp_format(mrep(2),3) ' (' num2str(100*CDFrep(2)) '%) ' ...
        smart_exp_format(mrep(3),3) ' (' num2str(100*CDFrep(3)) '%)'];
    nt=nt+1; titl{nt} = ['Magnitude statistics: ' ...
        smart_exp_format(mhi,3) ' (05%) ' ...
        smart_exp_format(mmd,3) ' (50%) ' ...
        smart_exp_format(mlo,3) ' (95%)'];
    nt=nt+1; titl{nt} = ' ';

    clf;
    ndx = X <= mmd;
    xx = (X(ndx)-mmd)/siglo/sqrt(2);
    Fmod = NaN(size(F));
    Fmod(ndx) = 1 - 0.5*(1+erf(xx));
    ndx = ~ndx;
    xx = (X(ndx)-mmd)/sighi/sqrt(2);
    Fmod(ndx) = 1 - 0.5*(1+erf(xx));
    plot([mlo mmd mhi],[0.95 0.50 0.05],'o', ...
        'MarkerFaceColor','r','MarkerEdgeColor','r');
    hold on;
    plot(X,Fmod,'--m','LineWidth',1);
    plot(X,F,'-k','LineWidth',1);
    hold off;
    
    set(gca, 'XDir', 'reverse');
    grid on;
    xlabl = ['Magnitude (Range = ' num2str(Rnormkm) ' km)'];
    ylabl = 'CDF';
    xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
    ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
    title(titl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
    set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);
    grid on;

    saveas(gcf,[fullfile(outpath,outfile) '_MagRNormCDF.png']);

    % Plot the range-norm magnitude
    clf;
    mrkr = 'o';
    lwid = 0.5;
    plot(dat.Phase,dat.MagNorm, ...
            'LineStyle','none', ...
            'LineWidth',lwid, ...
            'Marker',mrkr, ...
            'MarkerSize',msiz, ...
            'MarkerFaceColor','k', ...
            'MarkerEdgeColor','k');
    hold on;
    plot(xrng,[mrep(1) mrep(1)],'-.b','LineWidth',1);
    plot(xrng,[mrep(2) mrep(2)],'--b','LineWidth',1);
    plot(xrng,[mrep(3) mrep(3)],':b','LineWidth',1);
    plot(xrng,[mcut mcut],'-r','LineWidth',1);
    hold off;
    xlim(xrng);
    set(gca, 'YDir', 'reverse');
    yrng = plot_range([dat.MagNorm; mcut],0.05);
    ylim(yrng);
    xlabl ='Phase Angle (deg)';
    ylabl = ['Magnitude (Range = ' num2str(Rnormkm) ' km)'];
    xlabel(xlabl,'FontSize',xfsz,'FontWeight',xfwt);
    ylabel(ylabl,'FontSize',yfsz,'FontWeight',yfwt);
    
    % Put SATCON-1 cutoff in title
    nt=numel(titl)-1;
    titl{nt} = ['Fraction brighter than SATCON-1 limit: ' ...
                smart_exp_format(100*fcut,Ndig) '%'];
    title(titl,'FontSize',tfsz,'FontWeight',tfwt,'FontAngle',tfan);
    set(gca,'FontSize',afsz,'FontWeight',afwt,'LineWidth',alwd);

    saveas(gcf,[fullfile(outpath,outfile) '_MagPhase.png']);
    
end

% Save output structure
if UseMMTSaveFiles
    save([fullfile(outpath,outfile) '_AnalyzeMMTData.mat'],'out');
end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall    |  2024-Nov  | Initial version.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
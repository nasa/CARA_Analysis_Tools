function [Pc_all,Uc_all,Nc_tfc_all,Pc_att,Uc_att,Nc_tfc_att,Pc_0,Uc_0,Nc_0] = ...
    Pc_MC_Kep2body_parallel(Nsample_total, tfc_bins, r1, v1, P1, r2, v2, P2, ...
                            HBR, GM, motion_mode, conf_level, Nsample_batch, ...
                            plot_ca_dist_path)
% =========================================================================
%
% Function to estimate Pc with Monte Carlo sampling assuming Keplerian
% 2-body motion.
%
% =========================================================================
%
% INPUT:
%
%   Nsample_total       = Number of MC samples
%   tfc_bins            = Time of first contact bins
%   (r1,v1,P1)          = Primary state and cov
%   (r2,v2,P2)          = Secondary state and cov
%       NOTE: v1 and/or v2 = [] indicates equinoctial sampling, and r1/r2
%             hold the equinoctial states and P1/P2 hold the equinoctial
%             covariances
%   HBR                 = Combined pri+sec hard-body radii
%   GM                  = Central mass x grav constant
%   motion_mode         = Dynamical model: 'k2bpri', or 'k2balt'
%   conf_level          = Fractional confidence level for uncertainties, e.g., 
%                             0.6827 = chi2cdf(1,1) for 1-sigma interval
%                             0.95 for 95% interval
%                             0.999 for 99.9% interval
%                             chi2cdf(N^2,1) for N-sigma interval
%                             (optional, defaults to 1-sigma interval)
%   Nsample_batch       = Number of MC samples per parallel batch
%                         (optional, defaults to 1000 trials per processor)
%   plot_ca_dist_path   = Filepath at which to save Monte Carlo Output
%                         Plots (optional, does not print plots if no path
%                         specified)
%
% OUTPUT:
%
%   Pc_all              = Probability using all HBR time-of-first-contact
%                         (TFC) counts.
%   Uc_all              = Uncertainty on Pc_all.
%   Nc_tfc_all          = Number of all HBR TFC counts.
%   Pc_att              = Probability using only the first HBR TFC counts.
%                         This attenuated probability accounts for 
%                         destructive contacts.
%   Uc_att              = Uncertainty on Pc_att.
%   Nc_tfc_att          = Number of first HBR TFC counts.
%
% Examples/Validation Cases:
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required:     conj_bounds_Coppola.m
%                         refine_bounded_extrema.m
%                         Cart2Kep.m
%                         convert_cartesian_to_equinoctial.m
%                         jacobian_equinoctial_to_cartesian.m
%
% See also: None
%
% Revision History:
%   02/20/2018:           Added Accomodation of "0" state elements
%   07/08/2019:           Removed Rectilinear sampling schemes
% 
% Latest Revision:        07/08/2018
% =========================================================================

% Defaults

Nargin = nargin;

% Default to 1-sigma confidence interval
if (Nargin < 12) || isempty(conf_level) 
    conf_level = chi2cdf(1,1); 
end

% Default to 1000 samples per parallel batch
if (Nargin < 13) || isempty(Nsample_batch)
    Nsample_batch = 1e3; 
end

% Default to no distribution plotting
if (Nargin < 14)
    plot_ca_dist_path = [];
end
plot_ca_distribution = (~isempty(plot_ca_dist_path) && ~isnumeric(plot_ca_dist_path));

% Check if Velocity Covariance is populated
if size(P1,1)~=6
    error('Primary covariance matrix must be 6X6')
elseif sum(diag(P1)==0)~=0
    error('Primary covariance matrix must be fully populated')
elseif size(P2,1)~=6
    error('Secondary covariance matrix must be 6X6')
elseif sum(diag(P2)==0)~=0
    error('Secondary covariance matrix must be fully populated')
end

% Check if Statistics Toolbox Installed
use_statistics_toolbox = license('test','statistics_toolbox');

% Confidence level processing
if conf_level < 1 && conf_level > 0 && ~use_statistics_toolbox
    % Workaround for those without statistics toolbox to use sqrt(N) uncertainties
    binofit_uncertainties = false;
elseif conf_level < 1 && conf_level > 0 && use_statistics_toolbox
    % For 0 < ci < 1, use binofit to calculate asymmetric confidence 
    binofit_uncertainties = true;
else
    error('Confidence interval must be between 0 and 1');
end

% Number of sample processing
if (Nsample_batch < 1e2)
    % Do at least 100 samples per parallel batch    
    Nsample_batch = 1e2;
end

% Set up for binning
tbin = tfc_bins.xmd;
tmin = tfc_bins.x1;
tmax = tfc_bins.x2;

% Pad the bins so that the end bins span one half of a bin width
hdel = tfc_bins.del/2;
tbinpad = [tbin(1)-hdel tbin tbin(end)+hdel];

% Add (-Inf,Inf) to the bin ranges to enable histc to bin properly

binranges = [-Inf tfc_bins.xlo tfc_bins.xhi(end) Inf];

% Motion mode

motion_mode = lower(motion_mode);

switch motion_mode
    case {'k2bpri', 'k2b', 'kep2body', 'kep2bodypri', 'kep2body_primary'}
        linear_motion = false;
        use_primary_k2b_solver = true;
    case {'k2balt', 'kep2bodyalt', 'k2balt0', 'kep2body_alternate'}
        linear_motion = false;
        use_primary_k2b_solver = false;
    otherwise
        error('Invalid motion_mode parameter');
end        

% Calculate how many batches need to be calculated 

if (Nsample_total <= Nsample_batch)
    Nbatch = 1;
else
    fNbatch = Nsample_total/Nsample_batch;
    if rem(fNbatch,1) == 0
        Nbatch = fNbatch;
    else
        Nbatch = floor(fNbatch)+1;
    end
end

% Convert State Vectors to Equinoctial Frame if Required
% eqsamp1 = true; % Equinoctial sampling for primary
if ~isempty(v1)
    [KEP] = Cart2Kep([r1; v1]'/1000,'Mean','Rad'); % Get Keplerian elements to determine need for retrograde factor
    if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
        fr      = -1;
    else
        fr      = 1;
    end
    [~,n,af,ag,chi,psi,lM,~] = convert_cartesian_to_equinoctial(r1/1000,v1/1000,fr); % Equinoctial elements at TCA
    mu1    = [n af ag chi psi lM fr]';
    Jctoe1 = jacobian_equinoctial_to_cartesian(mu1,[r1; v1]'/1000,fr); % Jacobian going from equinoctial to cartesian
    Jctoe1 = inv(Jctoe1); % Jacobian going from cartesian to equinoctial

    P1 = Jctoe1 * (P1/1e6) * Jctoe1'; % Equinoctial covariance at TCA
    mu1    = mu1(1:6);
else
    mu1 = [r1(1:6); v1];
end

% Equinoctial sampling for secondary
if ~isempty(v2)
    [KEP] = Cart2Kep([r2; v2]'/1000,'Mean','Rad'); % Get Keplerian elements to determine need for retrograde factor
    if abs(KEP(3))>.95*pi && abs(KEP(3))<1.05*pi 
        fr      = -1;
    else
        fr      = 1;
    end
    [~,n,af,ag,chi,psi,lM,~] = convert_cartesian_to_equinoctial(r2/1000,v2/1000,fr); % Equinoctial elements at TCA
    mu2    = [n af ag chi psi lM fr]';
    Jctoe1 = jacobian_equinoctial_to_cartesian(mu2,[r2; v2]'/1000,fr); % Jacobian going from equinoctial to cartesian
    Jctoe1 = inv(Jctoe1); % Jacobian going from cartesian to equinoctial
    P2 = Jctoe1 * (P2/1e6) * Jctoe1'; % Equinoctial covariance at TCA
    mu2    = mu2(1:6);
else
    mu2 = [r2(1:6); v2];
end
    
% Set up for sampling using Eigen decomposition with nonnegative lambda
% clipping and include HBR uncertainties for both the primary and
% secondary covariances 

% Eigen decomp for Primary
[V1,D1] = eig(P1); 

L1 = diag(D1);
L1(L1 < 0) = 0;   % Make eigenvalues nonnegative

% Eigen decomp for secondary
[V2,D2] = eig(P2);

L2 = diag(D2);
L2(L2 < 0) = 0;   % Make eigenvalues nonnegative

% Define unused variables for Matlab's mvnrnd function
T1_mvnrnd = []; A1 = []; T2_mvnrnd = []; A2 = []; 
    
% Combined hard body radius squared
HBR2 = HBR^2;

% Set up for DCA minima and optimization
rbe_tolY = 1e-4*HBR2;
rbe_tolX = 1e-4*tfc_bins.del;
rbe_checks = false;
rbe_verbose = 0;

fzero_options = optimset('fzero');
fzero_options = optimset(fzero_options,'TolX',double(eps('single')));
fzero_options = optimset(fzero_options,'Display','notify');

% Buffers for distribution
if plot_ca_distribution
    TCAbuf = NaN(Nbatch,Nsample_batch);
    X1buf  = NaN(Nbatch,Nsample_batch,6); % Holds CA states
    X2buf  = NaN(Nbatch,Nsample_batch,6);
    Y1buf  = NaN(Nbatch,Nsample_batch,6); % Holds sampled states 
    Y2buf  = NaN(Nbatch,Nsample_batch,6);
end

% Initialize collision counters
Nbin = tfc_bins.N;
Nc_tfc_all = zeros(Nbatch,Nbin);
Nc_tfc_att = zeros(Nbatch,Nbin);
Nc_0 = zeros(Nbatch,1);

parfor nb=1:Nbatch
% for nb=1:Nbatch
    
    % Initialize "kep" structures, which will parameters for the primary
    % K2B solver.
    if use_primary_k2b_solver
        % Set "kep" structures to indicate primary K2B solver should
        % be used
        kep1.use_primary_k2b_solver = true;
        kep2.use_primary_k2b_solver = true;
    else
        % Set "kep" structures to indicate alternate K2B solver should
        % be used
        kep1.use_primary_k2b_solver = false;
        kep2.use_primary_k2b_solver = false;
    end    

    % Number of samples for this batch
    Nsample = min(Nsample_batch,Nsample_total-(nb-1)*Nsample_batch);
        
    % Sigma vectors for sampling primary and secondary PDFs
    S1 = sqrt(L1); S2 = sqrt(L2);

    % Sample using Eigen decomposition
    x1samp = repmat(mu1,[1 Nsample]) + [V1 * ...
        [S1(1)*randn(1,Nsample);  ...
         S1(2)*randn(1,Nsample);  ...
         S1(3)*randn(1,Nsample);  ...
         S1(4)*randn(1,Nsample);  ...
         S1(5)*randn(1,Nsample);  ...
         S1(6)*randn(1,Nsample)];  ...
         zeros(length(mu1)-6,Nsample)];
    x1samp = x1samp';

    x2samp = repmat(mu2,[1 Nsample]) + [V2 * ...
        [S2(1)*randn(1,Nsample);  ...
         S2(2)*randn(1,Nsample);  ...
         S2(3)*randn(1,Nsample);  ...
         S2(4)*randn(1,Nsample);  ...
         S2(5)*randn(1,Nsample);  ...
         S2(6)*randn(1,Nsample)];  ...
         zeros(length(mu1)-6,Nsample)];
    x2samp = x2samp';

    % Convert primary equinoctial states to ECI
    % Default to prograde elements if fr not specified
    if length(r1) == 7
        fr = r1(7);
    else
        fr = 1;
    end
    rvec = zeros(size(x1samp,1),3);
    vvec = rvec;
    for ii = 1:size(x1samp,1)
        [rvec(ii,:),vvec(ii,:),~,~,~,~,~,~,~] = ...
            convert_equinoctial_to_cartesian(x1samp(ii,1),x1samp(ii,2),x1samp(ii,3),x1samp(ii,4),x1samp(ii,5),x1samp(ii,6),0,fr);
    end
    x1samp = [rvec*1000 vvec*1000];
    % Convert secondary equinoctial states to ECI
    % Default to prograde elements if fr not specified
    if length(r2) == 7
        fr = r2(7);
    else
        fr = 1;
    end
    rvec = zeros(size(x2samp,1),3);
    vvec = rvec;
    for ii = 1:size(x2samp,1)
        [rvec(ii,:),vvec(ii,:),~,~,~,~,~,~,~] = ...
            convert_equinoctial_to_cartesian(x2samp(ii,1),x2samp(ii,2),x2samp(ii,3),x2samp(ii,4),x2samp(ii,5),x2samp(ii,6),0,fr);
    end
     x2samp = [rvec*1000 vvec*1000];
    
    % if plot_ca_distribution
        TCAbf = NaN(1,Nsample_batch);
        X1bf  = NaN(Nsample_batch,6);
        X2bf  = NaN(Nsample_batch,6);
        Y1bf  = NaN(Nsample_batch,6);
        Y2bf  = NaN(Nsample_batch,6);
    % end
    
    % Find the minimum of |r1(t)-r2(t)| for each sample
    for ns=1:Nsample
        % Currently sampled state vectors
        x1s = x1samp(ns,:)';
        x2s = x2samp(ns,:)';

        if use_primary_k2b_solver
            % Calculate the K2B solver parameters using t = 0.
            % This initializes the "kep1" and "kep2" structures.
            [kep1,~,x1s] = k2b_state_transition(x1s,0,GM,kep1);
            [kep2,~,x2s] = k2b_state_transition(x2s,0,GM,kep2);
        end

        % Anonymous function for distance squared, f(t) = D2 = |r2-r1|^2
        d2fun = @(t)kep_dist2(t,x1s,x2s,GM,kep1,kep2);
        
        % Calculate D2 at the padded bin-center time points
        d2binpad = d2fun(tbinpad);
        
        % Accumulate counts for overlaps at the initial time, tmin
        initial_overlap = (d2binpad(1) < HBR2);
        if initial_overlap
            Nc_0(nb) = Nc_0(nb) + 1;
        end

        % Find and refine all of the D2 minima within the range of the bins
        [tmma,d2mma,~,~,converged,~,tref,d2ref,~,~] = ...
            refine_bounded_extrema(d2fun,tbinpad,d2binpad,[],100,1, ...
                rbe_tolX,rbe_tolY,1,rbe_verbose,rbe_checks);

        if ~converged
            warning('Convergence not achieved in refine_bounded_extrema');
        end

        % Find if any minima imply incursions because D2 < HBR^2
        incur = d2mma < HBR2;
        Nincur = sum(incur);
        
        % Save for distribution
        if plot_ca_distribution
            [~,ndx] = min(d2mma);
            TCAb = tmma(ndx);
            TCAbf(ns) = TCAb; % TCA for first incursion
            X1b = kep_state(TCAbf(ns),x1s,GM,kep1);
            X1bf(ns,:) = X1b';
            X2b = kep_state(TCAbf(ns),x2s,GM,kep2);
            X2bf(ns,:) = X2b';
            Y1bf(ns,:) = x1s';
            Y2bf(ns,:) = x2s';
        end

        % Process incursions by finding the times of first contact, TFC,
        % defined by g(t) = D2 - HBR^2 = 0 and dg/dt < 0 (i.e., decreasing
        % roots).

        if (Nincur > 0)
            
            % Anonymous function for g(t) = D2 - HBR2
            ddfun = @(t)kep_ddif2(t,x1s,x2s,HBR2,GM,kep1,kep2);            

            % Find indices bounding all decreasing roots of g(t) = D2-HBR2
            nlo = find(diff(d2ref-HBR2 < 0) == 1);
            nhi = nlo+1;

            % Find all decreasing roots of g(t) = D2-HBR2, and increment
            % the incursion bin counters
            Nroot = numel(nlo);
            
            if (Nroot > 0)
                troot = NaN(1,Nroot);

                for nr=1:Nroot

                    % Bounding times for this root
                    tlo = tref(nlo(nr));
                    thi = tref(nhi(nr));

                    % Process if bounding times are within bin range
                    if (thi >= tmin) && (tlo <= tmax)

                        % Refine the bounded root
                        troot(nr) = fzero(ddfun,[tlo thi],fzero_options);

                    end

                end
                
                % Use histc to perform binning, neglecting the two ends
                % which are bounded by -Inf and +Inf, respectively
                if ~initial_overlap
                    
                    % Bin all of the roots (i.e., all of the incursions)
                    Nc_tfc = histc(troot,binranges);
                    Nc_tfc_all(nb,:) = Nc_tfc_all(nb,:) + Nc_tfc(2:end-2);

                    % Bin the first root for destructive incursions,
                    % which lead to attenuated counts
                    if (Nroot > 1)
                        Nc_tfc = histc(troot(1),binranges);
                    end
                    Nc_tfc_att(nb,:) = Nc_tfc_att(nb,:) + Nc_tfc(2:end-2);
                    
                end
                
            end

        end

    end
    
    if plot_ca_distribution
        TCAbuf(nb,:) = TCAbf; % TCA for first incursion
        X1buf(nb,:,:) = X1bf;
        X2buf(nb,:,:) = X2bf;
        Y1buf(nb,:,:) = Y1bf;
        Y2buf(nb,:,:) = Y2bf;
    end

    % Report info

    ts2 = current_timestring();    
    disp(['  End batch ' num2str(nb) ' of ' num2str(Nbatch) ' at ' ts2 ...
          ' (Nsample = ' num2str(Nsample) ')']);

end

% Count the overall number of collisions and calculate the probability
Nc_tfc_all = sum(Nc_tfc_all,1);
Nc_tfc_att = sum(Nc_tfc_att,1);

Nc_0 = sum(Nc_0,1);
if binofit_uncertainties
    [Pc_0,Uc_0] = binofit(Nc_0,Nsample_total,1-conf_level);
else
    Pc_0 = Nc_0/Nsample_total;    
    Uc_0 = sqrt(Nc_0)/Nsample_total;
end

% Nc_all = Nc_0 + sum(Nc_tfc_all);
Nc_all = sum(Nc_tfc_all);
if binofit_uncertainties
    [Pc_all,Uc_all] = binofit(Nc_all,Nsample_total,1-conf_level);
else
    Pc_all = Nc_all/Nsample_total;    
    Uc_all = sqrt(Nc_all)/Nsample_total;
end

% Nc_att = Nc_0 + sum(Nc_tfc_att);
Nc_att = sum(Nc_tfc_att);
if binofit_uncertainties
    [Pc_att,Uc_att] = binofit(Nc_att,Nsample_total,1-conf_level);
else
    Pc_att = Nc_att/Nsample_total;
    Uc_att = sqrt(Nc_att)/Nsample_total;
end

%% Generate Plots as Needed
% Plot the CA distributions

if plot_ca_distribution
    
    hfig = figure;
    
    % Calculate the distribution of close approaches
    Nca_max = Nbatch*Nsample_batch;
    
    Nca = 0;
    Xca = NaN(1,Nca_max);
    Yca = Xca;
    Tca = Xca;
    
    X1s = NaN(6,Nca_max);
    X2s = NaN(6,Nca_max);
    
    for nb=1:Nbatch
        ndx  = ~isnan(TCAbuf(nb,:));
        Nndx = sum(ndx);

        Tx = TCAbuf(nb,ndx);
        Dx = squeeze(X2buf(nb,ndx,:)-X1buf(nb,ndx,:));

        for nn=1:Nndx
            
            n = nn+Nca;
            
            zhat = Dx(nn,4:6)' / norm(Dx(nn,4:6));
            xhat = cross(zhat,[0; 0; 1]); xhat = xhat/norm(xhat);
            yhat = cross(zhat,xhat);
            Xca(n) = Dx(nn,1:3) * xhat;
            Yca(n) = Dx(nn,1:3) * yhat;
            Tca(n) = Tx(nn);
            
        end
        
        Nca1 = Nca+1;
        Nca2 = Nca+Nndx;
        
        X1s(:,Nca1:Nca2) = squeeze(X1buf(nb,ndx,:))';
        X2s(:,Nca1:Nca2) = squeeze(X2buf(nb,ndx,:))';
        
        Nca = Nca2;
        
    end
    
    if (Nca < Nca_max)
        Xca = Xca(1:Nca); Yca = Yca(1:Nca); Tca = Tca(1:Nca);
        X1s = X1s(:,1:Nca); X2s = X2s(:,1:Nca); 
    end
    
    Rca = sqrt(Xca.^2+Yca.^2);
    
    hits = Rca < HBR; Nhits = sum(hits);
    
    Ncirc = 3601;
    tcirc = linspace(0,2*pi,Ncirc);
    xcirc = HBR*cos(tcirc);
    ycirc = HBR*sin(tcirc);
    
    % Plot the distribution of close approaches, with a zoom
    
    Hzoom = 10;
    
    xmx0 = max([HBR max(abs(Xca)) max(abs(Yca))]);
    
    for npass=1:2
        
        nplt = 1+(npass-1)*3;
        subplot(2,3,nplt);
        
        if (npass == 1)
            % No zoom
            mrkr = '.';
            msiz = 1;
            xrng = plot_range([-xmx0 xmx0],0.05);
        else
            % Zoom in on HBR
            mrkr = '.';
            msiz = 3;
            xrng = Hzoom*[-HBR HBR];
        end
        
        % Scale
        
        if xrng(2) >= 1000
            Lscl = 1e3;
            Lstr = 'km';
        else
            Lscl = 1;
            Lstr = 'm';
        end
        
        % Plot all events in blue

        mcol = 'b';
        plot(Xca/Lscl,Yca/Lscl, ...
            'Marker',mrkr,'MarkerSize',msiz,'LineStyle','none', ...
            'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);
        
        hold on;

        % Plot hits in red

        mcol = 'r';
        plot(Xca(hits)/Lscl,Yca(hits)/Lscl, ...
            'Marker',mrkr,'MarkerSize',msiz,'LineStyle','none', ...
            'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);

        % Plot HBR boundary in green
        
        plot(xcirc/Lscl,ycirc/Lscl,'-g');
        
        % Lines through origin

        zrng = zeros(size(xrng));
        plot(xrng/Lscl,zrng/Lscl,'-k');
        plot(zrng/Lscl,xrng/Lscl,'-k');
        
        hold off;
        
        xlim(xrng/Lscl);
        ylim(xrng/Lscl);
        
        axis square;

        xlabel(['X_{CA} (' Lstr ')']);
        ylabel(['Y_{CA} (' Lstr ')']);
        
        drawnow;
        
    end
    
    % Create and plot the CDF(Rca)
    
    subplot(2,3,[2 3]);
    
    Rca_min = min(Rca);
    Rca_max = max(Rca);
    
    [FRcdf,Rcdf] = ecdf(Rca);
    
    if (Rca_min < HBR) && (Rca_max > HBR)
        RR = [Rcdf; HBR];
        FF = [FRcdf; interp1(Rcdf(2:end),FRcdf(2:end),HBR)];
        [Rcdf,ndx] = sort(RR);
        FRcdf = FF(ndx);
        clear RR FF;
    end
    
    xrng = 10.^plot_range(log10([HBR Rca_min Rca_max]),0.05);
    yrng = 10.^plot_range(log10(FRcdf(2:end)),0.05);
    
    FRcdf(1) = yrng(1);
    loglog(Rcdf,FRcdf,'-b');
    
    hold on;
    
    ndx = Rcdf <= HBR;
    plot(Rcdf(ndx),FRcdf(ndx),'-r');
    
    xlim(xrng);    
    ylim(yrng);
    
    xlabel('CA Distance (m)');
    ylabel('CDF');
    set(gca,'yaxislocation','right');        
    
    hold on;
    plot([HBR HBR],yrng,'-g')
    hold off;
    
    % Title
    titl = 'Equinoctial sampling at nominal TCA';
    titl = [titl ', HBR=' num2str(HBR) 'm'];
    title(titl,'FontWeight','bold','FontAngle','italic');
    
    drawnow;
    
    % Plot the distribution of CA distances in time
    
    subplot(2,3,[5 6]);
    xrng = plot_range(Tca,0.05);
    yrng = 10.^plot_range(log10([Rca_min Rca_max HBR]),0.05);
    
    mrkr = '.';
    msiz = 3;
    
    mcol = 'b';
    semilogy(Tca,Rca, ...
        'Marker',mrkr,'MarkerSize',msiz,'LineStyle','none', ...
        'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);

    hold on;
    
    mcol = 'r';
    plot(Tca(hits),Rca(hits), ...
        'Marker',mrkr,'MarkerSize',msiz,'LineStyle','none', ...
        'MarkerFaceColor',mcol,'MarkerEdgeColor',mcol);

    plot(xrng,[HBR HBR],'-g');
    
    hold off;

    xlim(xrng);
    ylim(yrng);
    
    xlabel('CA Relative Time (s)');
    ylabel('CA Distance');
    set(gca,'yaxislocation','right'); 
    
    if binofit_uncertainties
        [Pc_MC,Pc_MC_unc] = binofit(Nhits,Nca);
    else
        Pc_MC = Nhits/Nca;
        Pc_MC_unc = sqrt(Nhits)/Nca;
    end
    edel = max(abs(Pc_MC_unc-Pc_MC));
    [xstr , ~] = smart_error_format(Pc_MC,edel);
    [xstr1, ~] = smart_error_format(Pc_MC_unc(1),edel);
    [xstr2, ~] = smart_error_format(Pc_MC_unc(2),edel);
    
    Ncstr = smart_exp_format(Nhits,numel(num2str(Nhits,'%0.0f')));
    Nsstr = smart_exp_format(Nca,numel(num2str(Nca,'%0.0f')));
    if (Nca ~= str2double(Nsstr)); Nsstr = num2str(Nca); end
    titl = ['Nc=' Ncstr ', Ns=' Nsstr ', Pc=' xstr ...
            ' (95% ' xstr1 ' to ' xstr2 ')'];
    title(titl,'FontWeight','bold','FontAngle','italic');
    
    drawnow;
    
    fileroot = fullfile(plot_ca_dist_path,'CAdist');
    try
        saveas(hfig,[fileroot '.fig']);
    catch
        warning('Error Attempting to save .fig file, may be due to file size limitations')
    end
    saveas(hfig,[fileroot '.png']);
    
    % Plot the position and velocity distributions
    clf;
    
    mrkr = '.';
    msiz = 2;
    
    for nplt=1:4
        
        if nplt == 1
            xyz = X1s(1:3,:);
            titl = 'Pri.Pos.';
        elseif nplt == 2
            xyz = X1s(4:6,:);
            titl = 'Pri.Vel.';
        elseif nplt == 3
            xyz = X2s(1:3,:);
            titl = 'Sec.Pos.';
        elseif nplt == 4
            xyz = X2s(4:6,:);
            titl = 'Sec.Vel.';
        end
            
        subplot(2,2,nplt);
    
        mclr = 'b';    
        plot3(xyz(1,:),xyz(2,:),xyz(3,:), ...
            'Linestyle','none', ...
            'Marker',mrkr, ...
            'MarkerFaceColor',mclr, ...
            'MarkerEdgeColor',mclr, ...
            'MarkerSize',msiz);

        hold on;
        mclr = 'r';
        plot3(xyz(1,hits),xyz(2,hits),xyz(3,hits), ...
            'Linestyle','none', ...
            'Marker',mrkr, ...
            'MarkerFaceColor',mclr, ...
            'MarkerEdgeColor',mclr, ...
            'MarkerSize',msiz);
        hold off;
        
        % axis equal;
        axis tight;
        
        title(titl,'FontWeight','bold','FontAngle','italic');
        
        grid on;
        
        drawnow;
        
    end
    
    fileroot = fullfile(plot_ca_dist_path,'SAMPdist');
    % saveas(hfig,[fileroot '.fig']);
    saveas(hfig,[fileroot '.png']);
    
    delete(hfig); clear hfig
    
    save([fileroot '.mat'],'X1s','X2s','hits');
  
end

return;
end

% =========================================================================

function dist2 = kep_dist2(t,x10,x20,GM,kep1,kep2)

% Calculate the distance-squared between two objects

% Calculate the states
x1 = kep_state(t,x10,GM,kep1);
x2 = kep_state(t,x20,GM,kep2);

% Calculate the squares of the separation distances 
dist2 = sum((x1(1:3,:)-x2(1:3,:)).^2,1);

return;
end

% =========================================================================

function ddif2 = kep_ddif2(t,x10,x20,HBR2,GM,kep1,kep2)
% Kep. 2-body separation distance squared minus HBR^2
ddif2 = kep_dist2(t,x10,x20,GM,kep1,kep2)-HBR2;

return;
end

% =========================================================================

function x1 = kep_state(t,x10,GM,kep1)

% Calculate the Kep.2-body or linear-motion state for time t
if isempty(kep1)
    
    % Linear motion
    st = size(t);
    r1 = repmat(x10(1:3),st);
    v1 = repmat(x10(4:6),st);
    tt = repmat(t,[3 1]);
    x1 = [r1 + tt .* v1; v1];
    
else
    % Kepler 2-body motion
    if kep1.use_primary_k2b_solver
        % Primary K2B state calculator, "K2Bpri" or just "K2B"
        % disp('Using primary K2B state calculator');
        [~,~,x1] = k2b_state_transition(x10,t,GM,kep1);
        
    else
        % Alternate K2B state calculator
        % disp('Using alternate K2B state calculator');
        Nt = numel(t);
        x1 = zeros(6,Nt);
        
        % Using new versions of Cart2Kep and Kep2Cart, optimized for this
        % task, "K2Balt"
        [aa,ee,ii,OO,ww,MM] = Cart2Kep_MeanAnom_MKS(x10(1:3),x10(4:6));
        
        nn = sqrt(GM/aa^3);
        Mt = MM + nn.*t;
        
        for j=1:Nt
            [x1(1:3,j),x1(4:6,j)] = Kep2Cart_MeanAnom_MKS(aa,ee,ii,OO,ww,Mt(j));
        end
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
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 08/07/2019 | Modification From Original Analysis Code
%                               Base to Simplify Code

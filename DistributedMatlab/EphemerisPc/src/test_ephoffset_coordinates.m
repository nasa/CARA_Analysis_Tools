function test_ephoffset_coordinates(eph)
% test_ephoffset_coordinates - Test ephemeris offset concepts using the
%                              input ephemeris table.
%
% Syntax: test_ephoffset_coordinates(eph);
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
%    eph        -   Ephemeris structure with fields:
%                     .T             - Ephemeris times      [1xN]
%                     .X             - State vectors        [6xN]
%                     .P             - Covariance matrices  [6x6xN]
%                     .N             - Number of points
%                     .SecPerTimeUnit - Seconds per time
%                                      unit (optional,
%                                      default = 86400)
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Initialize random numbers
rng default;

% If not given in eph. table, assume eph. time units of days
if ~isfield(eph,'SecPerTimeUnit'); eph.SecPerTimeUnit = []; end
if isempty(eph.SecPerTimeUnit); eph.SecPerTimeUnit = 86400; end

% Median time between eph points
dTmedian = median(diff(eph.T));

% Do some testing using the primary and secondary eph tables
rmag = NaN(1,eph.N);
for i=1:eph.N
    rmag(i) = norm(eph.X(1:3,i));
end

% Function for radial distance
rmagfun = @(tt) interpRadialDist(tt,eph);

% Refine the minima and maxima
Ttol = [1e-2*dTmedian NaN]; % Absolute and relative time tolerance
Rtol = [1e-1 NaN]; % Absolute and relative radial distance tolerance (km)
extrema_types = 3; % Refine both minima and maxima
endpoints = false; % Refine at endpoints
rbeverbose = false;
rbecheck = false;
[trmin,rmin,trmax,rmax,rbeconv,rbenbisect] = ...
        refine_bounded_extrema(rmagfun,eph.T,rmag, ...
                               [],100,extrema_types,Ttol,Rtol, ...
                               endpoints,rbeverbose,rbecheck);
if ~rbeconv
    error(['radial distance minima search failed to converge in ' ...
        num2str(rbenbisect) ' bisections']);
end

figure(1); clf;
plot(NaN,NaN); hold on;
xrng = plot_range(eph.T,0.05);
yrng = plot_range(rmag,0.05,0.5);
xlim(xrng); ylim(yrng);
plot(trmin,rmin,'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(trmax,rmax,'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(eph.T,rmag,'.k');
hold off;

% Combine minima and maxima into one list of times, and add some
% intervening times as well

tlist = unique([trmin trmax]);
Nlist = numel(tlist);
tadd = [];
uadd = (1:3)/4;
for n=1:Nlist-1
    ta = tlist(n);
    tb = tlist(n+1);
    tadd = [tadd ta+uadd*(tb-ta)]; %#ok<AGROW>
end
tlist = unique([eph.T(1) eph.T(end) tlist tadd]);
Nlist = numel(tlist);

% Calculate cov realism stats
Nsample = 1e2;
Q0 = 1/(12*Nsample);
chicdfmod = (2*(1:Nsample)-1)./(2*Nsample);
I6x6 = eye(6,6); I3x3 = eye(3,3);
twopi = 2*pi;

% Extract eph info
QE = NaN(size(tlist)); QX = QE; QY = QE;
QXr = QE; QYr = QE; QXv = QE; QYv = QE;
SXr = NaN(3,Nlist); SXv = SXr; SYr = SXr; SYv = SXr;
for i = 1:Nlist
    
    % Extract time
    t = tlist(i);
    
    % Interpolate mean cartesian state and cov
    [r,v,PX,na,~,~,~,nb,~,~,~] = interpStateCov(t,eph);
    X = [r v]'; invPX = PX\I6x6;
    invPXr = PX(1:3,1:3)\I3x3; invPXv = PX(4:6,4:6)\I3x3;
    
    % Convert cartesian state/cov to equinoctial state/cov
    [~,nmean,afmean,agmean,chimean,psimean,lMmean,~] = ...
        convert_cartesian_to_equinoctial(r,v);
    lMmean = mod(lMmean,twopi);    
    E = [nmean,afmean,agmean,chimean,psimean,lMmean]';

    % Calculate Jacobian matrix and inverse
    J = jacobian_equinoctial_to_cartesian(E,X); % dX/dE
    K = J\I6x6; % dE/dX

    % Calculate equinoctial covariance, enforcing symmetry to correct
    % round-off error effects
    PE = cov_make_symmetric( K * PX * K' ); invPE = PE\I6x6;

    % Sample eq states
    Es = mvnrnd(E,PE,Nsample);
    % Maha distances squared
    dEs = Es' - repmat(E,[1 Nsample]);
    M2 = sum(dEs.*(invPE*dEs),1);
    % Sorted M2 values
    M2 = sort(M2);
    % CvM statistic 
    Qcdf = chicdfmod-chi2cdf(M2,6);
    QE(i) = Q0 + sum(Qcdf(:).^2);
    
    % Sample cart states
    Xs = NaN(size(Es));
    for s=1:Nsample
        [rs,vs] = convert_equinoctial_to_cartesian( ...
            Es(s,1),Es(s,2),Es(s,3),Es(s,4),Es(s,5),Es(s,6),0);
        Xs(s,:) = [rs' vs'];
    end
    
    % Maha distances squared
    dXs = Xs' - repmat(X,[1 Nsample]);
    M2 = sum(dXs.*(invPX*dXs),1);
    % Sorted M2 values
    M2 = sort(M2);
    % CvM statistic 
    Qcdf = chicdfmod-chi2cdf(M2,6);
    QX(i) = Q0 + sum(Qcdf(:).^2);
    
    % Maha distances squared
    dXs = Xs(:,1:3)' - repmat(X(1:3),[1 Nsample]);
    M2 = sum(dXs.*(invPXr*dXs),1);
    % Sorted M2 values
    M2 = sort(M2);
    % CvM statistic 
    Qcdf = chicdfmod-chi2cdf(M2,3);
    QXr(i) = Q0 + sum(Qcdf(:).^2);
    % Sigma values
    [~,D] = eig(PX(1:3,1:3)); SXr(:,i) = sqrt(sort(diag(D)));
    
    % Maha distances squared
    dXs = Xs(:,4:6)' - repmat(X(4:6),[1 Nsample]);
    M2 = sum(dXs.*(invPXv*dXs),1);
    % Sorted M2 values
    M2 = sort(M2);
    % CvM statistic 
    Qcdf = chicdfmod-chi2cdf(M2,3);
    QXv(i) = Q0 + sum(Qcdf(:).^2);
    % Sigma values
    [~,D] = eig(PX(4:6,4:6)); SXv(:,i) = sqrt(sort(diag(D)));
    
    % Sample ephemeris offset states
    Ys = NaN(size(Es));
    for s=1:Nsample
        [YY,C2E] = convert_cartesian_to_ephoffset(t,Xs(s,:)',X,eph,1);
        if any(isnan(YY))
            break;
        end
        Ys(s,:) = YY';
    end
    
    % If all samples are defined, then process further
    if ~any(isnan(Ys(:)))
        
        [Y,C2E] = convert_cartesian_to_ephoffset(t,X,X,eph,3);
        % JYapp = jacobian_cartesian_to_ephoffset(t,X,X,eph); % dY/dX
        M3   = [C2E.Xhat_tau  ; C2E.Yhat_tau  ; C2E.Zhat_tau];
        JY = [M3 zeros(3,3); zeros(3,3) M3];
        PY = cov_make_symmetric( JY * PX * JY' );
        invPY = PY\I6x6; invPYr = PY(1:3,1:3)\I3x3; invPYv = PY(4:6,4:6)\I3x3;
        
        % Maha distances squared
        dYs = Ys' - repmat(Y,[1 Nsample]);
        M2 = sum(dYs.*(invPY*dYs),1);
        % Sorted M2 values
        M2 = sort(M2);
        % CvM statistic for state 
        Qcdf = chicdfmod-chi2cdf(M2,6);
        
        % Maha distances squared
        QY(i) = Q0 + sum(Qcdf(:).^2);
        dYs = Ys(:,1:3)' - repmat(Y(1:3),[1 Nsample]);
        M2 = sum(dYs.*(invPYr*dYs),1);
        % Sorted M2 values
        M2 = sort(M2);
        % CvM statistic for position
        Qcdf = chicdfmod-chi2cdf(M2,3);
        QYr(i) = Q0 + sum(Qcdf(:).^2);
        % Sigma values
        [~,D] = eig(PY(1:3,1:3)); SYr(:,i) = sqrt(sort(diag(D)));
        
        % Maha distances squared
        dYs = Ys(:,4:6)' - repmat(Y(4:6),[1 Nsample]);
        M2 = sum(dYs.*(invPYv*dYs),1);
        % Sorted M2 values
        M2 = sort(M2);
        % CvM statistic 
        Qcdf = chicdfmod-chi2cdf(M2,3);
        QYv(i) = Q0 + sum(Qcdf(:).^2);
        % Sigma values
        [~,D] = eig(PY(4:6,4:6)); SYv(:,i) = sqrt(sort(diag(D)));
        
    end

end

% Index values of rmag minima and maxima
imin = ismember(tlist,trmin);
imax = ismember(tlist,trmax);

% Plot limits
Qcmb = [QX QE QY]; Qmin = min(Qcmb); Qmax = max(Qcmb);
yrng = 10.^plot_range(log10([Qmin Qmax]),0.05);
xrng = plot_range(tlist,0.05);

% Filled area for realistic covariances (CCW from bottom left)
Qlimit = 1.162;
xfill = [xrng(1) xrng(2) xrng(2) xrng(1) xrng(1)];
yfill = [yrng(1) yrng(1) Qlimit  Qlimit  yrng(1)];
cfill = 0.7*[1 1 1];

% Make plots of Q stats for full 6 dimensional states
figure(2); clf;
subplot(3,1,1);
semilogy(NaN,NaN); hold on;
fill(xfill,yfill,cfill);
plot(tlist(imin),QE(imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(tlist(imax),QE(imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(tlist,QE,'.k');
hold off;
xlim(xrng); ylim(yrng);
% xlabel('Time (days)');
ylabel('Q (equinoctial)');

subplot(3,1,2);
semilogy(NaN,NaN); hold on;
fill(xfill,yfill,cfill);
plot(tlist(imin),QX(imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(tlist(imax),QX(imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(tlist,QX,'.k');
hold off;
xlim(xrng); ylim(yrng);
% xlabel('Time (days)');
ylabel('Q (cartesian)');

subplot(3,1,3);
semilogy(NaN,NaN); hold on;
fill(xfill,yfill,cfill);
plot(tlist(imin),QY(imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(tlist(imax),QY(imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(tlist,QY,'.k');
hold off;
xlim(xrng); ylim(yrng);
xlabel('Time (days)');
ylabel('Q (eph.offset)');
drawnow;

% Make plots of Q stats for 3 dimensional position states
figure(3); clf;

Qcmb = [QXr QXv QYr QYv]; Qmin = min(Qcmb); Qmax = max(Qcmb);
yrng = 10.^plot_range(log10([Qmin Qmax]),0.05);

subplot(2,2,1);
semilogy(NaN,NaN); hold on;
fill(xfill,yfill,cfill);
plot(tlist(imin),QXr(imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(tlist(imax),QXr(imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(tlist,QXr,'.k');
hold off;
xlim(xrng); ylim(yrng);
% xlabel('Time (days)');
ylabel('Q (cartesian position)');

subplot(2,2,3);
semilogy(NaN,NaN); hold on;
fill(xfill,yfill,cfill);
plot(tlist(imin),QYr(imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(tlist(imax),QYr(imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(tlist,QYr,'.k');
hold off;
xlim(xrng); ylim(yrng);
xlabel('Time (days)');
ylabel('Q (eph.offset position)');

subplot(2,2,2);
semilogy(NaN,NaN); hold on;
for i=1:3
    plot(tlist(imin),SXr(i,imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(tlist(imax),SXr(i,imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
    plot(tlist,SXr(i,:),'.k');
end
hold off;
xlim(xrng); % ylim(yrng);
% xlabel('Time (days)');
ylabel('\sigma (cartesian position)');

subplot(2,2,4);
semilogy(NaN,NaN); hold on;
for i=1:3
    plot(tlist(imin),SYr(i,imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(tlist(imax),SYr(i,imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
    plot(tlist,SYr(i,:),'.k');
end
hold off;
xlim(xrng); % ylim(yrng);
xlabel('Time (days)');
ylabel('\sigma (eph.offset position)');
drawnow;

% Make plots of Q stats for 3 dimensional position states
figure(4); clf;

subplot(2,2,1);
semilogy(NaN,NaN); hold on;
fill(xfill,yfill,cfill);
plot(tlist(imin),QXv(imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(tlist(imax),QXv(imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(tlist,QXv,'.k');
hold off;
xlim(xrng); ylim(yrng);
% xlabel('Time (days)');
ylabel('Q (cartesian velocity)');

subplot(2,2,3);
semilogy(NaN,NaN); hold on;
fill(xfill,yfill,cfill);
plot(tlist(imin),QYv(imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
plot(tlist(imax),QYv(imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
plot(tlist,QYv,'.k');
hold off;
xlim(xrng); ylim(yrng);
xlabel('Time (days)');
ylabel('Q (eph.offset velocity)');

subplot(2,2,2);
semilogy(NaN,NaN); hold on;
for i=1:3
    plot(tlist(imin),SXv(i,imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(tlist(imax),SXv(i,imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
    plot(tlist,SXv(i,:),'.k');
end
hold off;
xlim(xrng); % ylim(yrng);
% xlabel('Time (days)');
ylabel('\sigma (cartesian velocity)');

subplot(2,2,4);
semilogy(NaN,NaN); hold on;
for i=1:3
    plot(tlist(imin),SYv(i,imin),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(tlist(imax),SYv(i,imax),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
    plot(tlist,SYv(i,:),'.k');
end
hold off;
xlim(xrng); % ylim(yrng);
xlabel('Time (days)');
ylabel('\sigma (eph.offset velocity)');
drawnow;

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
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
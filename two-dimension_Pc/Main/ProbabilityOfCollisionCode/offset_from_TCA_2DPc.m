function [out] = offset_from_TCA_2DPc(X1TCA,C1TCA,X2TCA,C2TCA,HBR,params)
% offset_from_TCA_2DPc_ephem - Generate offset-from-TCA 2D-Pc estimates for
%                              a conjunction, with offset times spanning
%                              the conjunction's short-term encounter
%                              validity interval.
%
% Syntax: [out] = offset_from_TCA_2DPc_ephem(X1TCA,C1TCA,X2TCA,C2TCA,HBR,params)
%
% As described by Hall (2019), a conjunction that satisfies rigorously the
% 2D-Pc method assumptions of linear-motion and constant pos. covariance
% has small variations in offset-from-TCA 2DPc values over the short-term
% encounter validity interval (STEVI)
%      TCA - dtstv <= time <= TCA + dtstv
% with the STEVI half-width time dtstv given by Coppola (2012b) eq (17).
% This function calculates offset-from-TCA 2D-Pc values over the STEVI for
% both raw, unremediated equinoctial-state TCA covariances and remediated
% equinoctial-state TCA covariances.
%
% Inputs:
%
%   X1TCA       = Primary cartesian TCA state (m) [6x1]
%   C1TCA       = Primary cartesian TCA covariance (m^2 & (m/s)^2) [6x6]
%   X2TCA       = Secondary cartesian TCA state (m) [6x1]
%   C2TCA       = Secondary cartesian TCA covariance (m^2 & (m/s)^2) [6x6]
%   HBR         = Combined hard-body radius (m)
%   params      = Auxilliary input parameters (optional):
%     (NOTE: See "Initializations and Defaults" section below for
%      definitions and details of the auxilliary parameters.)
%
% Outputs:
%
%   out.any_remediation = Flag indicating that primary and/or secondary
%                         TCA equinoctial state covariances required NPD
%                         remediation.
%   out.Unrem = Strtucture holding results for unremediated covariances
%       Unrem.tau0 = Conjunction duration begin time relative to TCA
%       Unrem.tau1 = Conjunction duration end   time relative to TCA
%       Unrem.dtstv = Short-term encounter validity interval half-width
%       Unrem.toff = Initial offset-from-TCA time array relative to TCA
%       Unrem.Pc2Doff = Offset-from-TCA 2DPc values at the toff times
%       Unrem.tref = Refined offset-from-TCA time array relative to TCA
%       Unrem.Pc2Dref = Offset-from-TCA 2D-Pc values at the tref times
%       Unrem.Pc2Dmnma = Number of offset-from-TCA Pc2D minima in the
%                        STEVI interval -dtstv < t-TCA < dtstv
%       Unrem.Pc2Dmxma = Number of offset-from-TCA Pc2D maxima in the
%                        STEVI interval -dtstv < t-TCA < dtstv
%       Unrem.imnma = Indices for the offset-from-TCA Pc2D minima for the
%                     tref and Pc2Dref arrays
%       Unrem.imnma = Indices for the offset-from-TCA Pc2D maxima for the
%                     tref and Pc2Dref arrays
%   out.Remed = Strtucture holding quantities for remediated covariances
%       (with contents the of same form as described above for out.Unrem)
%   out.PriEq = Structure holding primary TCA states and covariances
%       PriEq.E = Equinoctial state [n,af,ag,chi,psi,lM]'
%       PriEq.Jetoc = Jacobian matrix, dX/dE
%       PriEq.Jctoe = Inverse of Jacobian matrix
%       PriEq.Praw = Raw TCA equinoctial state covariance matrix
%       PriEq.Lraw = Raw TCA eq.state covariance eigenvalues
%       PriEq.Vraw = Raw TCA eq.state covariance eigenvector matrix
%       PriEq.Lrem = Remediated TCA eq.stae covariance eigenvalues
%       PriEq.Prem = Remediated TCA equinoctial state covariance matrix
%       PriEq.Prem = Determinant of Prem
%       PriEq.Pinv = Inverse of Prem
%       PriEq.Pinv = Inverse of Prem
%       PriEq.Crem = Remediated TCA cartesian state (pos,vel) covariance
%       PriEq.Pposdefstat = PD status of Praw
%           -1 => Praw is NPD
%            0 => Praw is PSD
%           +1 => Praw is PD
%       PriEq.Pclipstat = Clipping status of Praw
%           true  => Min eigenvalue(s) of Lraw clipped to create Lrem
%           false => Min eigenvalue(s) of Lraw not clipped, so Lrem = Lraw
%   out.SecEq = Structure holding primary TCA states and covariances
%       (with contents of the same form as described above for out.PriEq)
%
% References:
%
%    V.T. Coppola (2012a) "Including Velocity Uncertianty in the
%    Probability of Collision Between Space Objects" AAS 12-247.
%
%    V.T. Coppola (2012b) "Evaluating the Short Encounter Assumption
%    of the Probability of Collision Formula" AAS 12-248.
%
%    D.T.Hall (2019) "Implementation Recommendations and Usage Boundaries
%    for the Two-Dimensional Probability of Collision Calculation"
%    AAS 19-632.
%
% Example/Validation Cases:
%
% Other m-files required:
%   conj_bounds_Coppola.m
%   convert_cartesian_to_equinoctial.m
%   CovRemEigValClip.m
%   jacobian_equinoctial_to_cartesian.m
%   refine_bounded_extrema.m
%
% Subfunctions:
%   EqTCAStCov
%   CalcSingleOffset
%   CalcMultiOffset
%
% MAT-files required: None
%
% See also: None
%
% September 2019; Last revision: 2019-SEP-18
%
% ----------------- BEGIN CODE -----------------

% Initializations and Defaults

Nargin = nargin;

if Nargin < 6; params = []; end

% 2D-Pc algorithm - defaults to PcElrod.m 2D-Pc function
if ~isfield(params,'Pc2DAlg') || isempty(params.Pc2DAlg)
    params.Pc2DAlg = 'ELR'; % Alternatives include 'FE92' and 'AA00'
end

% Min number of offset-from-TCA ephemeris points
if ~isfield(params,'Noffmin') || isempty(params.Noffmin)
    params.Noffmin = 101;
end
% Make Noffmin an odd number, so that it has a point at the center of the
% STEVI, which coincides with the TCA itself
if mod(params.Noffmin,2) == 0
    params.Noffmin = params.Noffmin+1;
end

% Max number of offset-from-TCA ephemeris points
if ~isfield(params,'Noffmax') || isempty(params.Noffmax)
    params.Noffmax = 1001; % Max number of offset time points
end
% Make Noffmax an odd number
if mod(params.Noffmax,2) == 0
    params.Noffmax = params.Noffmax-1;
end

% Precision factor for conjunction durations
if ~isfield(params,'gamma') || isempty(params.gamma)
    params.gamma = 1e-16; % Coppola AAS 12-248 conj. duration parameter
end

% X tolerance for refining extrema in the offset-from-TCA 2D-Pc curve
if ~isfield(params,'tolXfact') || isempty(params.tolXfact)
    params.tolXfact = 1e-6;
end

% Y tolerance for refining extrema in the offset-from-TCA 2D-Pc curve
if ~isfield(params,'tolYfact') || isempty(params.tolYfact)
    params.tolYfact = 1e-3;
end

% Verbosity
if ~isfield(params,'verbose') || isempty(params.verbose)
    params.verbose = false;
end

% Other initializations

twopi = 2*pi;

% Parameters required to search for minima and maxima in the
% offset-from-TCA Pc curve

Nbisectmax = 100;    % Maximum number of bisections in extrema search
extrema_types = 3;   % Search for both minima and maxima in search
endpoints = false;   % Search for only interior extrema - not edge extrema
check_inputs = true; % Check inputs on extrema search

% Calculate the primary and secondary equinoctial TCA states and
% covariances, and remediate the covariances

out.PriEq = EqTCAStCov(X1TCA,C1TCA); E1TCA = out.PriEq.E;
out.SecEq = EqTCAStCov(X2TCA,C2TCA); E2TCA = out.SecEq.E; 

% NPD status for TCA 6x6 equinoctial state covariance matrices

if params.verbose
    disp(['Equinoctial TCA covariance PD status for pri. and sec.: ' ...
        num2str(out.PriEq.Pquality) ' ' num2str(out.SecEq.Pquality)]);
    disp(['Equinoctial TCA cov. clipping status for pri. and sec.: ' ...
        num2str(out.PriEq.Pclipstat) ' ' num2str(out.SecEq.Pclipstat)]);
end

% Flag if any TCA equinoctial covariance remediation occured

out.any_remediation = out.PriEq.Pclipstat | out.SecEq.Pclipstat;

% Calculate the orbital periods using TCA equinoctial states, 
% and the mininumum half-period

period1 = twopi/out.PriEq.E(1);
period2 = twopi/out.SecEq.E(1);
hperiod = min(period1,period2)/2;
    
if params.verbose
    disp([' Orbital period for pri & sec, and min(period)/2 = ' ...
        num2str(period1) ' ' ...
        num2str(period2) ' ' ...
        num2str(hperiod) ' s']);
end

% Calculate unremediated Coppola (2012) conjunction encounter duration
% using the relative TCA state and covariance

rci = X2TCA(1:3)-X1TCA(1:3);         % Rel position
vci = X2TCA(4:6)-X1TCA(4:6);         % Rel velocity
Pci = C2TCA(1:6,1:6)+C1TCA(1:6,1:6); % Rel covariance
[out.Unrem.tau0,out.Unrem.tau1] = ...
    conj_bounds_Coppola(params.gamma,HBR,rci,vci,Pci);

if params.verbose
    disp(['Conjunction bounds, midpoint and duration: ' ...
        num2str(out.Unrem.tau0) ' ' num2str(out.Unrem.tau1) ' ' ...
        num2str(0.5*(out.Unrem.tau1+out.Unrem.tau0)) ' ' ...
        num2str(out.Unrem.tau1-out.Unrem.tau0)]);
end

if out.any_remediation
    % Eigenvalue-clipping remediation means possible change in duration
    Pci = out.SecEq.Crem(1:6,1:6)+out.PriEq.Crem(1:6,1:6);
    [out.Remed.tau0,out.Remed.tau1] = ...
        conj_bounds_Coppola(params.gamma,HBR,rci,vci,Pci);
    if params.verbose
        disp(['Remediated conjunction bounds, midpoint and duration: ' ...
            num2str(out.Remed.tau0) ' ' num2str(out.Remed.tau1) ' ' ...
            num2str(0.5*(out.Remed.tau1+out.Remed.tau0)) ' ' ...
            num2str(out.Remed.tau1-out.Remed.tau0)]);
    end
else
    % No Eigenvalue-clipping remediation means no difference in duration
    out.Remed.tau0 = out.Unrem.tau0; out.Remed.tau1 = out.Unrem.tau1;
end

% Calculate the offset-from-TCA 2D-Pc ephemeris for both the unremediated
% and remediated TCA equinoctial covariances, using one or two passes as
% necessary

if out.any_remediation
    Npass = 2; % Two passes needed for remediated and unremediated cases
else
    Npass = 1; % Only one pass needed because no remediation occurred
end

for npass=1:Npass
    
    % Conjunction encounter duration
    
    if npass == 1
        % Use unremediated conj. duration and TCA equinoctial covariances
        tau0  = out.Unrem.tau0;
        tau1  = out.Unrem.tau1;
        P1TCA = out.PriEq.Praw;
        P2TCA = out.SecEq.Praw;
    else
        % Use remediated conj. duration and TCA equinoctial covariances
        tau0  = out.Remed.tau0;
        tau1  = out.Remed.tau1;
        P1TCA = out.PriEq.Prem;
        P2TCA = out.SecEq.Prem;
    end
    
    % Short-term encounter validity interval half-width
    % from Coppola (2012b) eq (17); also see Hall (2019).

    dtstv = max([tau1-tau0 abs(tau0) abs(tau1)]);
    
    % Generate the limits of the offset times to explore. Specifically,
    % span Coppola's entire short-term encounter validity interval (STEVI):
    %
    %        TCA-dtstv <= time <= TCA+dtstv
    % or
    %           -dtstv <= toff <= +dtstv
    
    toffmin = -dtstv; 
    toffmax =  dtstv;
    
    % Define offset time array
    
    toffdif = toffmax-toffmin;
    
    if toffdif <= hperiod
        % STEVI less than or equal to min orbital period => use the min
        % number of offset time points for Noff
        Noff = params.Noffmin;
    else
        % STEVI greater than to min orbital period => increase
        % proportionately the number of offset time points up to the
        % maximum limit, and ensure Noff is an odd number
        Noff0 = 2*round((params.Noffmin-1)*hperiod/toffdif/2)+1;
        Noff = min(params.Noffmax,max(params.Noffmin,Noff0));
    end
    
    toff = linspace(toffmin,toffmax,Noff);

    % Generate the 2D-Pc values over all of the offset times
    
    Pc2Doff = CalcMultiOffset(toff,E1TCA,P1TCA,E2TCA,P2TCA,HBR, ...
                              params.Pc2DAlg,params.verbose);
    
    % Extract the at-TCA 2D-Pc value, which is the midpoint of the curve
                    
    Nmid = (Noff-1)/2+1;
    Pc2DTCA = Pc2Doff(Nmid);

    % Refine the offset-Pc2D extrema, using the function
    % refine_bounded_extrema.m

    % Tolerances for extrema refinement convergence
    tolX = params.tolXfact * toffdif;
    tolY = params.tolYfact * max(Pc2Doff);
    
    % Anonymous function for extrema refinement
    Pc2Dfun = @(t)CalcMultiOffset(t,E1TCA,P1TCA,E2TCA,P2TCA,HBR, ...
                                  params.Pc2DAlg,params.verbose);

    % Perform extrema refinement
    [xmnma,ymnma,xmxma,ymxma,cnv,nbs,tref,Pc2Dref,imnma,imxma] =    ...
        refine_bounded_extrema(Pc2Dfun,toff,Pc2Doff,[],Nbisectmax,  ...
            extrema_types,tolX,tolY,endpoints,                      ...
            params.verbose,check_inputs); %#ok<ASGLU>

    % % Debugging code to plot results of extrema search
    % figure;
    % subplot(2,1,1);
    % plot(toff,Pc2Doff,'-k');
    % hold on;
    % plot(tref,ymnma,'*r');
    % plot(xmxma,ymxma,'*r');
    % hold off;
    % xrng = plot_range(toff,0.05);
    % xlim(xrng);
    % subplot(2,1,2);    
    % semilogy(toff,Pc2Doff,'-k');
    % hold on;
    % plot(xmnma,ymnma,'*r');
    % plot(xmxma,ymxma,'*r');
    % hold off;
    % xlim(xrng);
    % keyboard;

    % Save the required output quantities
    
    if npass == 1
        % Save output variables calculated with the unremediated covariance
        out.Unrem.dtstv    = dtstv;
        out.Unrem.toff     = toff;
        out.Unrem.Pc2Doff  = Pc2Doff;
        out.Unrem.Pc2DTCA  = Pc2DTCA;
        out.Unrem.tref     = tref;
        out.Unrem.Pc2Dref  = Pc2Dref;
        out.Unrem.Pc2Dmnma = numel(imnma);
        out.Unrem.Pc2Dmxma = numel(imxma);
        out.Unrem.imnma    = imnma;
        out.Unrem.imxma    = imxma;
        % If no cov. remediation occured, then the remediated output is the
        % same as the unremediated output (and 2nd pass is not performed)
        if ~out.any_remediation
            out.Remed.dtstv    = dtstv;
            out.Remed.toff     = toff;
            out.Remed.Pc2Doff  = Pc2Doff;
            out.Remed.Pc2DTCA  = Pc2DTCA;
            out.Remed.tref     = tref;
            out.Remed.Pc2Dref  = Pc2Dref;
            out.Remed.Pc2Dmnma = numel(imnma);
            out.Remed.Pc2Dmxma = numel(imxma);
            out.Remed.imnma    = imnma;
            out.Remed.imxma    = imxma;
        end
    else
        % Save output variables calculated with the remediated covariance
        % during the 2nd pass
        out.Remed.dtstv    = dtstv;
        out.Remed.toff     = toff;
        out.Remed.Pc2Doff  = Pc2Doff;
        out.Remed.Pc2DTCA  = Pc2DTCA;
        out.Remed.tref     = tref;
        out.Remed.Pc2Dref  = Pc2Dref;
        out.Remed.Pc2Dmnma = numel(imnma);
        out.Remed.Pc2Dmxma = numel(imxma);
        out.Remed.imnma    = imnma;
        out.Remed.imxma    = imxma;
    end
    
end

return
end

% =========================================================================

function [out] = EqTCAStCov(XTCA,CTCA)

% Calculate the TCA equinoctial state and covariance, and remediate the
% covariance using eigenvalue clipping

X = XTCA/1e3; % Cartesian state (km & km/s) [6x1]
PX = CTCA/1e6; % Covariance (km & km/s) [6x6]
[~,n,af,ag,chi,psi,lM,~] = ...
    convert_cartesian_to_equinoctial(X(1:3),X(4:6));
out.E = [n;af;ag;chi;psi;lM]; % TCA equinoctial state [6x1]
out.Jetoc = jacobian_equinoctial_to_cartesian(out.E,X);
out.Jctoe = out.Jetoc\eye(size(out.Jetoc));
out.Praw = out.Jctoe * PX * out.Jctoe'; % TCA equinoctial covariance
[out.Lrem,out.Lraw,out.Vraw,out.Pposdefstat,out.Pclipstat, ...
    out.Pdet,out.Pinv,out.Prem] = CovRemEigValClip(out.Praw,0);
PXrem = out.Jetoc * out.Prem * out.Jetoc';
out.Crem = PXrem * 1e6;  % Remediated cartesian state cov (m & m/s)^2    

if (min(diag(out.Praw)) < 0)
    out.Pquality = -2;
else
    out.Pquality = out.Pposdefstat;
end

return
end

% =========================================================================

function [PcOffset] = CalcSingleOffset(X1,C1,X2,C2,HBR,Pc2DAlg)

% Estimate an offset value using states and covariances that correspond to 
% a time that is potentially offset from the conjunction's nominal TCA

% Input primary and secondary cartesian states, which are not necessarily
% at the nominal TCA

r1 = X1(1:3); v1 = X1(4:6);
r2 = X2(1:3); v2 = X2(4:6);

% Calculate the "effective" TCA, assuming linear relative motion starting
% from the input primary and secondary cartesian states

r = r2-r1; v = v2-v1;
TCAeff = -(r'*v)/(v'*v);

% Primary and secondary states at the effective TCA

r1eff = r1+TCAeff*v1;
r2eff = r2+TCAeff*v2;

% Calculate 2DPc values using the effective TCA primary+secondary states
% along with the ephemeris covariances

switch Pc2DAlg
    
    case 'FE92' % Foster and Estes (1992)
        
        PcOffset = Pc2D_Foster( ...
            r1eff',v1',C1(1:3,1:3),    ...
            r2eff',v2',C2(1:3,1:3),    ...
            HBR,1e-9,'circle');
        
    case 'AA00' % Akella and Alfriend (2000)
        
        PcOffset = Pc_AkellaAlfriend( ...
            r1eff,v1,C1(1:3,1:3),     ...
            r2eff,v2,C2(1:3,1:3),     ...
            HBR);
        
    case 'ELR' % Elrod vectorized version of FE92 algorithm
        
        PcOffset = PcElrod(         ...
            r1eff',v1',C1(1:3,1:3), ...
            r2eff',v2',C2(1:3,1:3), ...
            HBR,64);
        
    otherwise
        
        error('Invalid Pc2DAlg parameter');

end

return
end

% =========================================================================

function [Pc2Doff] = CalcMultiOffset(toff,E1TCA,P1TCA,E2TCA,P2TCA,HBR,Pc2DAlg,verbose)

% Estimate offset values using states and covariances that correspond to 
% times that are potentially offset from the conjunction's nominal TCA

% Initialize the offset equinoctial states

E1 = E1TCA; E2 = E2TCA;

% Initialize the equinoctial 2-body motion STMs

S1 = eye(6,6);  S2 = S1;

% Allocate arrays for the offset-from-TCA 2DPc-related quantities

Soff = size(toff);
Noff = prod(Soff);

Pc2Doff = NaN(Soff);

% Loop over offset-from-TCA ephemeris points

for n=1:Noff

    % Propagate the states and covariances to the offset time using
    % the VA15 equinoctial element algorithm for Kepler 2-body
    % (K2B) motion 

    E1(6) = E1TCA(6)+E1TCA(1)*toff(n); % MeanLong = MeanLong0 + t*MeanMotion)
    S1(6,1) = toff(n); % STM(6,1) = d(MeanLong)/d(MeanMotion)
    P1 = S1 * P1TCA * S1'; % Propagate equinoctial covariance from TCA
    [r1,v1] = convert_equinoctial_to_cartesian( ...
        E1(1),E1(2),E1(3),E1(4),E1(5),E1(6),0);
    X1 = [r1;v1]; % Cartesian state (km,km/s)
    Jetoc1 = jacobian_equinoctial_to_cartesian(E1,X1); % Transformation matrix from E to X
    C1 = Jetoc1 * P1 * Jetoc1'; % Convert equinoctial cov into cartesian cov

    E2(6) = E2TCA(6)+E2TCA(1)*toff(n);
    S2(6,1) = toff(n);
    P2 = S2 * P2TCA * S2';
    [r2,v2] = convert_equinoctial_to_cartesian( ...
        E2(1),E2(2),E2(3),E2(4),E2(5),E2(6),0);
    X2 = [r2;v2];
    Jetoc2 = jacobian_equinoctial_to_cartesian(E2,X2);
    C2 = Jetoc2 * P2 * Jetoc2';

    % Convert states and covs from (km,km/s) to to (m,m/s) units

    X1 = X1*1e3; C1 = C1*1e6;
    X2 = X2*1e3; C2 = C2*1e6;
    
    % Calculate the 2DPc using the offset-from-TCA states and
    % covariances

    Pc2Doff(n) = CalcSingleOffset(X1,C1,X2,C2,HBR,Pc2DAlg);

    if verbose
        disp([num2str(toff(n),'%+0.3g') '  ' ...
            num2str(Pc2Doff(n),'%0.6g') ]);
    end

end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D.Hall         | 2019-SEP-17 | Initial Development
%
function [PcMax,PcMin,Nc,out] = EphPc(eph1,eph2,HBR,params)

% EphPc - Estimate collision probability bounds, expected number of
%         collisions, and collision rate using the 3D-Nc method and
%         ephemeris interpolation.
%
% Syntax: [PcMax, PcMin, Nc, out] = EphPc(eph1, eph2, HBR, params);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Estimates the collision probability bounds (PcMax and PcMin), the
%   statistically expected number of collisions (Nc) and the collision
%   rate using the 3D-Nc method and ephemeris interpolation.
%   Specifically, this function interpolates states & covariances to a
%   time t, used as initial conditions (i.e., t = t0 = t10 = t20).
%   Interpolation is done using the closest and bracketing eph points,
%   yielding 4 combinations that are used to calculate weighted means
%   and variances by blending.
%
% =========================================================================
%
% Input:
%
%    eph1       -   Primary ephemeris structure with fields:
%                     .T             - Ephemeris times      [1xN]
%                     .X             - State vectors        [6xN]
%                     .P             - Covariance matrices  [6x6xN]
%                     .N             - Number of points
%                     .SecPerTimeUnit - Seconds per time unit
%
%    eph2       -   Secondary ephemeris structure (same fields as eph1)
%
%    HBR        -   Combined hard body radius (km)
%
%    params     -   Input parameter structure (see set_default_param
%                   calls in code for full list of parameters and
%                   defaults)
%
% =========================================================================
%
% Output:
%
%   PcMax       -   Upper limit collision probability
%
%   PcMin       -   Lower limit collision probability
%
%   Nc          -   Statistically expected number of collisions
%
%   out         -   Output structure with fields:
%
%     Primary Outputs:
%     
%     .Nc          - Statistically expected number of collisions
%                    accumulated over the risk assessment interval
%     .Nccum       - Cumulative Nc (aligned with time grid .T)
%     .Ncdot       - Nc rate (aligned with .T)
%     .Ncseg       - Nc for each encounter segment
%     .PcMax       - Max collision probability accumulated over the
%                    risk assessment interval.
%     .PcMin       - Min collision probability accumulated over the
%                    risk assessment interval
%     .PcumMax     - Max cumulative Pc (aligned with .T)
%     .PcumMin     - Min cumulative Pc (aligned with .T)
%
%     Additional Outputs:
%     .cseg        - Index in .T of each .Ncdot (Nc rate) peak 
%     .eph1        - Struct containing primary ephemeris data
%     .eph2        - Struct containing secondary ephemeris data
%     .EPHBEG      - Ephemeris start MATLAB date number
%     .epos        - End index of each encounter segment determined by the 
%                    first zero in .Ncdot bracketing the segment on the
%                    right side
%     .eseg        - End index of each encounter segment in .T based on the
%                    .Ncdot minima bracketing the segment on the right side
%     .HBR         - Hard body radius (km)
%     .ipos        - Start index of each encounter segment determined by 
%                    the first zero in .Ncdot bracketing the segment on the
%                    left side
%     .iseg        - Start index of each encounter segment in .T based on 
%                    the .Ncdot minima bracketing the segment on the left
%                    side
%     .MDeff       - Effective Mahalanobis distance squared (aligned with 
%                    .T)
%     .MScut0      - Initial min collision-sphere Mahalanobis distance 
%                    cutoff used
%     .Nccmn       - Cumulative min Nc (aligned with .T)
%     .Nccmx       - Cumulative max Nc (aligned with .T)
%     .Ncdmn       - Min Nc rate (aligned with .T)
%     .Ncdmx       - Max Nc rate (aligned with .T)
%     .Ncmn        - Min Nc accumulated over the risk assessment interval
%     .Ncmx        - Max Nc accumulated over the risk assessment interval
%     .Ncsmn       - Min Nc for each encounter segment
%     .Ncsmx       - Max Nc for each encounter segment
%     .NT          - Number of time samples in .T
%     .Nseg        - Number of encounter segments
%     .ref         - Reference structure:
%                      .MD2      Mean-State Mahalanobis distance squared at
%                                the center of the collisioin sphere
%                                (aligned with .ref.T)
%                      .MD2max   MD2 at maxima (aligned with .ref.TMD2max)
%                      .MD2min   MD2 at minima (aligned with .ref.TMD2min)
%                      .MS2      Approximated minimum Mahalanobis distance 
%                                squared on the collision sphere 
%                                (aligned with .ref.T)
%                      .MS2min   MS2 values at minima (aligned with 
%                                .ref.TMS2min)
%                      .NT       Number of time points 
%                      .RD2      Squared relative distance (aligned with 
%                                .ref.T)
%                      .RD2min   RD2 values at minima (aligned with 
%                                .ref.TRD2min)
%                      .T        Time grid
%                      .TMD2max  Times of MD2 maxima 
%                      .TMD2min  Times of MD2 minima 
%                      .TMS2min  Times of MS2 minima 
%                      .TRD2min  Times of RD2 minima 
%     .Scdot       - Small HBR limit collision rate approximation (aligned
%                    with .T)
%     .tau_clipped - Flag for clipping risk assessment interval to 
%                    the intersection of ephemerides
%     .taua        - Risk assessment interval start time (days)
%     .taua0       - Risk assessment interval initial start time (days)
%     .taub        - Risk assessment interval end time (days)
%     .taub0       - Risk assessment interval initial end time (days)
%     .T           - Time grid in days from risk assessment interval start            
%     .Uc          - Nc interpolation uncertainty 
%     .Uccum       - Cumulative Nc interpolation uncertainty
%     .Ucseg       - Nc interpolation uncertainty for each segment
%     .UseEqEls    - Flag to use equinoctial elements
%     .Ucdot       - Uncertainty of Nc rate
%     .Vcdot       - Weighted uncertainty of Nc rate
%
% =========================================================================
%
% References:
%
%   Hall, D. T. "Ephemeris-Based Satellite Collision Rates and
%   Probabilities," Journal of Spacecraft and Rockets, Vol. 62, No. 4,
%   pp. 1152-1169, 2025.
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

% Set default parameters
params = set_default_param(params, 'MScut' , 1e3);
params = set_default_param(params, 'taua'  , -Inf);
params = set_default_param(params, 'taub'  ,  Inf);
params = set_default_param(params, 'Nctottol' , 5.0e-3);
params = set_default_param(params, 'Ncsegtol' , 5.0e-2);
params = set_default_param(params, 'Tcsegtol' , 2.5e-2);
params = set_default_param(params, 'Nctnytol' , 1.0e-1);
params = set_default_param(params, 'Tctnytol' , 1.0e-1);
params = set_default_param(params, 'RetainOneSegment' , true);
params = set_default_param(params, 'bracketing' , true);
params = set_default_param(params, 'Fclip' , 1e-4);
params = set_default_param(params, 'initial_refined_times' , 1);
params = set_default_param(params, 'RefineRelDistMinima' , false);
params = set_default_param(params, 'POPmaxiter' , 200);
params = set_default_param(params, 'StepFracMin' , 0.1);
params = set_default_param(params, 'NcFracTiny' , 1e-10);
params = set_default_param(params, 'NcFracNegligible' , 1e-30);
params = set_default_param(params, 'use_Lebedev' , true);
params = set_default_param(params, 'deg_Lebedev' , 5810);
params = set_default_param(params, 'wgt_Lebedev' , []);
params = set_default_param(params, 'vec_Lebedev' , []);
params = set_default_param(params, 'EPHBEG' , []);
params = set_default_param(params, 'tauwarning' , false);
params = set_default_param(params, 'SaveInitialSearch' , true);
params = set_default_param(params, 'MakeOutFile' , true);
params = set_default_param(params, 'outputRoot' , '');
params = set_default_param(params, 'outputPath' , pwd);
params = set_default_param(params, 'UseEqEls' , true);
params = set_default_param(params, 'CreateCCSDSOEMFiles' , false);
params = set_default_param(params, 'MSZeroNcdot' , params.MScut);
params = set_default_param(params, 'verbose' , true);

% Extract risk assessment interval bounds
taua = params.taua; taub = params.taub;

% Initialize output
PcMax      = NaN;   % Upper limit Pc calculated using Eph3DNc method
PcMin      = NaN;   % Lower limit Pc calculated using Eph3DNc method
Nc         = NaN;   % Total Nc       calculated using Eph3DNc method
out.PcMax  = PcMax;
out.PcMin  = PcMin;
out.Nc     = Nc;
out.Uc     = NaN;   % Nc interp. uncertainty estimated using Eph3DNc method
out.HBR    = HBR;   % HBR value used in Eph3DNc method (km)
out.taua0  = taua;  % Risk assessment initial start time 
out.taub0  = taub;   % Risk assessment initial start time 
out.MScut0 = params.MScut;  % Init. min. collision-sphere Maha.dist. cutoff
out.EPHBEG = params.EPHBEG; % Begin Matlab date number for ephemeris
out.UseEqEls = params.UseEqEls; % Flag to use equinoctial elements

% Initialize misc parameters
verbose = params.verbose;
debug_plotting = 0;

% Probability cutoffs for defining statistically significant Ncdot
% maxima
probcut_iter = [0.999 0.99 0.05];
probcut_last = [0.999 0.99 0.95];

% % Initialize testing mode for ephemeris-offset coordinates
% ephoffset_testing = false;

% Use persistent variables to prevent repetitive recalculation
% of Lebedev quadrature vectors and weights
persistent sph_Lebedev deg_Lebedev
if params.use_Lebedev && isempty(params.vec_Lebedev)
    if isempty(sph_Lebedev)
        % Calculate Lebedev quadrature quantities for the first time
        if verbose
            disp(['Calculating ' num2str(params.deg_Lebedev) ...
                  ' Lebedev quadrature points']);
        end
        sph_Lebedev = getLebedevSphere(params.deg_Lebedev);
        deg_Lebedev = params.deg_Lebedev;
    elseif deg_Lebedev ~= params.deg_Lebedev
        % Recalculate Lebedev quadrature quantities
        sph_Lebedev = getLebedevSphere(params.deg_Lebedev);
        deg_Lebedev = params.deg_Lebedev;
        % Issue warning, if appropriate
        if Lebedev_warning && ~params.suppress_Lebedev_warning
            warning(['Lebedev quadrature points recalculated; ' ...
                     'needless repetition can slow execution.']);
        end
    end
    % Copy Lebedev vector and weights from persistent variable
    params.vec_Lebedev = [sph_Lebedev.x sph_Lebedev.y sph_Lebedev.z]';
    params.wgt_Lebedev = sph_Lebedev.w;
end

% Create the CCSDS OEM files, if required
if params.CreateCCSDSOEMFiles > 0
    if isempty(params.EPHBEG)
        warning('EPHBEG parameter must be specified to create CCSDS OEM files');
    else
        ForceCreation = params.CreateCCSDSOEMFiles > 1;
        Success = CreateOEMFile(1, eph1, params.EPHBEG, params.outputPath, ForceCreation);
        if ~Success
            warning('Success not achieved in creating OEM file for primary');
        end
        Success = CreateOEMFile(2, eph2, params.EPHBEG, params.outputPath, ForceCreation);
        if ~Success
            warning('Success not achieved in creating OEM file for secondary');
        end
    end
end

% Process null risk assessment intervals
if (taua >= taub)
    PcMax = 0; out.PcMax = 0;
    PcMin = 0; out.PcMin = 0;
    Nc    = 0; out.Nc    = 0; out.Uc = 0;
    return;
end

% Ensure ephemeris tables have at least five points for interpolations
if eph1.N < 5 || eph2.N < 5
    error('Both ephemeris table need at least five time points');
end

% Ensure the ephmemeris tables are monotonically increasing in time
moninc1 = all(diff(eph1.T) > 0);
moninc2 = all(diff(eph2.T) > 0);
if ~moninc1 || ~moninc2
    error('Ephemeris table(s) must be monotonically increasing');
end

% Ensure the ephemeris tables overlap in time
if eph1.T(1) >= eph2.T(eph2.N) || eph1.T(eph1.N) <= eph2.T(1)
    error('Ephemeris tables do not overlap in time');
end

% Ensure both ephemeris tables span the risk assessment period
newa = max([eph1.T(1) eph2.T(1) taua]);
newb = min([eph1.T(eph1.N) eph2.T(eph2.N) taub]);
if newa ~= taua || newb ~= taub
    out.tau_clipped = true;
    wrnstr = 'Risk assessment interval clipped to intersection of ephemerides';
    if verbose
        disp(wrnstr);
    end
    if params.tauwarning 
        warning(wrnstr);
    end
    taua = newa;
    taub = newb;
else
    out.tau_clipped = false;
end
out.taua = taua; out.taub = taub;

% Report ephemeris parameters
if verbose
        disp('Primary ephemeris parameters;');
        disp([' Span (days): ' num2str(eph1.T(end)-eph1.T(1))]);
        disp([' Number of points: ' num2str(eph1.N)]);
        difTimes = diff(eph1.T)*eph1.SecPerTimeUnit;
        disp([' Min/med/max time steps (s): ' ...
              num2str(min(difTimes)) '  ' ...
              num2str(median(difTimes)) '  ' ...
              num2str(max(difTimes))]);
        disp('Secondary ephemeris parameters;');
        disp([' Span (days): ' num2str(eph2.T(end)-eph2.T(1))]);
        disp([' Number of points: ' num2str(eph2.N)]);
        difTimes = diff(eph2.T)*eph2.SecPerTimeUnit;
        disp([' Min/med/max time steps (s): ' ...
              num2str(min(difTimes)) '  ' ...
              num2str(median(difTimes)) '  ' ...
              num2str(max(difTimes))]);
        num1in2 = sum(ismember(eph1.T,eph2.T));
        num2in1 = sum(ismember(eph2.T,eph1.T));
        disp(['Number of pri. eph. points in sec. eph: ' num2str(num1in2)])
        disp(['Number of sec. eph. points in pri. eph: ' num2str(num2in1)])
end

% Eigenvalue clipping parameter (km^2)
Lclip = (params.Fclip*HBR)^2;

% Min step-size fraction to allow when combining ephemeris times
StepFracMin = params.StepFracMin;

% Tolerances for adaptively-refined trapezoidal integration
% (using iterative time-step bisection method)
Nctottol = params.Nctottol;
Ncsegtol = params.Ncsegtol; Tcsegtol = params.Tcsegtol;
Nctnytol = params.Nctnytol; Tctnytol = params.Tctnytol;

% Check for a previously saved mat file containing refined MD curves
if params.SaveInitialSearch
    SaveFile = [fullfile(params.outputPath,params.outputRoot) '_InitSearch.mat'];
    if exist(SaveFile,'file')
        outsav = out;
        if verbose
            disp(' ');
            disp('Loading save mat file');
            load(SaveFile,'out');
        end
        if ~isfield(out,'taua0')  || out.taua0  ~= outsav.taua0   || ...
           ~isfield(out,'taub0')  || out.taub0  ~= outsav.taub0   || ...
           ~isfield(out,'MScut0') || out.MScut0 ~= outsav.MScut0
            out = outsav;
            if verbose
                disp('Obsolete reference curve mat file');
            end
        end
    end
else
    warning('SaveInitialSearch parameter set to false');
end

% Check if out structure already has refined Maha. and relative distance
% reference curves. If not, generate a new reference curve mat file
if ~isfield(out,'ref') || ~isfield(out,'Ncdot') || ...
   (params.RefineRelDistMinima && ~isfield(out.ref,'RD2min'))

    % Find limits of pri and sec ephemeris times spanning (taua,taub)
    ndx = find((taua <= eph1.T) & (eph1.T <= taub));
    if isempty(ndx)
        [~,na] = min(abs(taua-eph1.T));
        [~,nb] = min(abs(taub-eph1.T));
        ndx = unique([na nb]);
    end
    n1a = ndx(1);
    n1b = ndx(end);
    ndx = find((taua <= eph2.T) & (eph2.T <= taub));
    if isempty(ndx)
        [~,na] = min(abs(taua-eph2.T));
        [~,nb] = min(abs(taub-eph2.T));
        ndx = unique([na nb]);
    end
    n2a = ndx(1);
    n2b = ndx(end);

    % Define the initial times for the refined MD curves
    if params.initial_refined_times == 1
        % Use primary eph times
        T = eph1.T(n1a:n1b);
    elseif params.initial_refined_times == 2
        % Use secondary eph times
        T = eph2.T(n2a:n2b);
    else
        % Combine pri & sec times, but eliminate points that are too close
        T = combine_eph_times(eph1.T(n1a:n1b),eph2.T(n2a:n2b), ...
            params.StepFracMin,verbose);
    end

    % Add the risk assessment interval end-point times
    T = combine_eph_times([taua taub],T,StepFracMin,verbose);

    % Ensure there are enough initial time points
    NTmin = 21;
    dT = [diff(T) 0];
    dTmax = (T(end)-T(1))/NTmin;
    ndx = dT > dTmax;
    while any(ndx)
        Tnew = sort([T T(ndx)+dT(ndx)/2]);
        T = Tnew;
        dT = [diff(T) 0];
        ndx = dT > dTmax;
    end

    % Calculate the separation and Mahalanobis distances
    if verbose
        disp(' ');
        disp('Calculating initial separation and Mahalanobis distance curves');
    end
    [MS2,MD2,~] = interp_MS2_MD2_PD2(T,eph1,eph2,HBR,Lclip);

    if debug_plotting > 0
        figure(1); clf;
        subplot(4,1,1);
        plot(T,MS2,'d-k','MarkerSize',4);
        hold on;
        plot(T,MD2,'o-c','MarkerSize',4);    
        hold off;
        ylabel('M2');
    end 
    
    % Refine the minima of the MS2 curve
    if verbose
        disp(' ');
        disp('Finding the extrema of Mahalanobis distance curves');
    end

    % Function for the Maha. distance at the center of the collision sphere
    MD2fun = @(TT) interp_MS2_MD2_PD2(TT,eph1,eph2,0,Lclip);
    
    % Function for approximated minimum Maha. dist. on the collision sphere
    MS2fun = @(TT) interp_MS2_MD2_PD2(TT,eph1,eph2,HBR,Lclip);
    
    % Parameters for bisection search to find the minimal
    Ttol = [NaN NaN];     % Impose no time tolerance
    MD2tol = [1e-3 1e-3]; % Impose both absolute and relative MD2 tolerance
    extrema_types = 1; % Refine only minima (not maxima)
    endpoints = true; % Refine at endpoints
    rbeverbose = verbose;
    rbecheck = true;
    
    % Find the minima of the MS2 curve
    % (i.e., the approximated minimum Maha. dist. on the collision sphere)
    [out.ref.TMS2min,out.ref.MS2min,~,~, ...
        rbeconv,rbenbisect,Ttmp,MS2tmp,~,~] = ...
            refine_bounded_extrema(MS2fun,T,MS2, ...
                                   [],100,extrema_types,Ttol,MD2tol, ...
                                   endpoints,rbeverbose,rbecheck);
    if ~rbeconv
        error(['MS2 minima search failed to converge in ' ...
            num2str(rbenbisect) ' bisections']);
    end

    if debug_plotting > 0
        subplot(4,1,2);
        plot(T,MS2,'xm','MarkerSize',4);
        hold on;
        plot(Ttmp,MS2tmp,'-k');
        plot(out.ref.TMS2min,out.ref.MS2min,'vr');
        hold off;
        ylabel('MS2');
    end
    
    % Find the minima of the MD2 curve
    % (i.e., the Maha. distance at the center of the collision sphere)
    MD2tmp = MD2fun(Ttmp);
    [out.ref.TMD2min,out.ref.MD2min,~,~, ...
        rbeconv,rbenbisect,out.ref.T,out.ref.MD2,~,~] = ...
            refine_bounded_extrema(MD2fun,Ttmp,MD2tmp, ...
                                   [],100,extrema_types,Ttol,MD2tol, ...
                                   endpoints,rbeverbose,rbecheck);
    if ~rbeconv
        error(['MD2 minima search failed to converge in ' ...
            num2str(rbenbisect) ' bisections']);
    end
    out.ref.MS2 = MS2fun(out.ref.T);

    if debug_plotting > 0
        subplot(4,1,3);
        plot(T,MD2,'xm','MarkerSize',4);
        hold on;
        plot(out.ref.T,out.ref.MD2,'-c');
        plot(out.ref.TMD2min,out.ref.MD2min,'vr');
        hold off;
        ylabel('MD2');
    end    

    % Find the maxima of the MD2 curve
    % (i.e., the Maha. distance at the center of the collision sphere)
    MD2tol = [NaN 1e-3]; % Impose only relative MD2 tolerance
    extrema_types = 2; % Refine only maxima (not minima)
    [~,~,out.ref.TMD2max,out.ref.MD2max,...
        rbeconv,rbenbisect,~,~,~,~] = ...
            refine_bounded_extrema(MD2fun,out.ref.T,out.ref.MD2, ...
                                   [],100,extrema_types,Ttol,MD2tol, ...
                                   endpoints,rbeverbose,rbecheck);
    if ~rbeconv
        error(['MD2 maxima search failed to converge in ' ...
            num2str(rbenbisect) ' bisections']);
    end

    if debug_plotting > 0
        subplot(4,1,4);
        plot(out.ref.T,out.ref.MD2,'-c');
        hold on;
        plot(out.ref.T,out.ref.MS2,'-k');
        plot(out.ref.TMD2min,out.ref.MD2min,'vr');
        plot(out.ref.TMS2min,out.ref.MS2min,'vm');
        plot(out.ref.TMD2max,out.ref.MD2max,'^b');
        % plot(out.ref.TMS2max,out.ref.MS2max,'^c');
        hold off;
        ylabel('M2');
        drawnow;
    end    

    % Add the maxima to the refined curve
    [out.ref.T,ndx] = sort([out.ref.T out.ref.TMD2max]);
    out.ref.MD2 = [out.ref.MD2 out.ref.MD2max];
    out.ref.MD2 = out.ref.MD2(ndx);
    out.ref.MS2 = [out.ref.MS2 MS2fun(out.ref.TMD2max)];
    out.ref.MS2 = out.ref.MS2(ndx);
    [out.ref.T,ndx] = unique(out.ref.T);
    out.ref.MD2 = out.ref.MD2(ndx);
    out.ref.MS2 = out.ref.MS2(ndx);
    
    % Refine the relative distance minima, only if required
    if params.RefineRelDistMinima
        
        if verbose
            disp(' ');
            disp('Refining the minima of the relative distance curve');
        end
        
        PD2fun = @(TT) interp_PD2(TT,eph1,eph2);
        Ttol = [NaN NaN];   % Impose no time tolerance
        RDtol = [1e-6 NaN]; % Impose both absolute and relative tolerance
        extrema_types = 1;  % Refine minima
        endpoints = false;  % Refine at endpoints
        rbeverbose = false;
        rbecheck = true;
        
        [out.ref.TRD2min,out.ref.RD2min,~,~, ...
            rbeconv,rbenbisect,new_ref_T,out.ref.RD2,~,~] = ...
                refine_bounded_extrema(PD2fun,out.ref.T,PD2fun(out.ref.T), ...
                                       [],100,extrema_types,Ttol,RDtol, ...
                                       endpoints,rbeverbose,rbecheck);
        if ~rbeconv
            error(['PD2 minima search failed to converge in ' ...
                num2str(rbenbisect) ' bisections']);
        end
        
        % Calculate new MD2 and MS2 points for output reference curves
        ndx = ismember(new_ref_T,out.ref.T); idx = ~ndx;
        tmp = out.ref.MD2; % Old points
        out.ref.MD2 = NaN(size(new_ref_T));
        out.ref.MD2(ndx) = tmp;
        out.ref.MD2(idx) = MD2fun(new_ref_T(idx));
        tmp = out.ref.MS2;
        out.ref.MS2 = NaN(size(new_ref_T));
        out.ref.MS2(ndx) = tmp; % Old points
        out.ref.MS2(idx) = MS2fun(new_ref_T(idx));
        out.ref.T = new_ref_T;
    
    end

    % Define number of points in refined MD curve
    out.ref.NT = numel(out.ref.MS2);

    if debug_plotting > 0
        figure(2); clf;
        plot(out.ref.T,out.ref.MD2,'-c');
        hold on;
        plot(out.ref.T,out.ref.MS2,':k');
        plot(out.ref.TMS2min,out.ref.MS2min,'vr');
        plot(out.ref.TMD2max,out.ref.MD2max,'^b');
        hold off;
    end

    % Determine if any of the refined MS2 curve points are below the
    % specified threshold.  If not, then set Nc = 0 and return.
    % If so, then add the points with MS2 <= MS2cut to the curves.
    % If all of the points are below the cutoff, do nothing.
    MS2cut = params.MScut^2;
    ndx = find(out.ref.MS2 <= MS2cut);
    Nndx = numel(ndx);

    if (Nndx == 0)

        % All MS2 curve points are above cutoff, so by inference Nc = 0
        Nc    = 0; out.Nc    = 0; out.Uc = 0;
        PcMax = 0; out.PcMax = 0;
        PcMin = 0; out.PcMin = 0;
        out.T = [taua taub]; out.Ncdot = [0 0]; out.Ucdot = [0 0];
        return;

    elseif (Nndx == out.ref.NT)

        % All MS2 curve points are below the cutoff 
        out.T = T;

    else

        % Some but not all MS2 curve points are below the cutoff, 
        % so add the points with MS2 = MS2cut to the refined curve
        MS2above = out.ref.MS2-MS2cut;
        
        % Refine the time segments defined by the MS2 <= MS2cut condition
        % First, find the curve indices that bound the segment begin points
        ilo = find(diff(MS2above < 0) == 1);
        ihi = ilo+1;

        % Find the curve indices that bound the segment end points
        elo = find(diff(MS2above > 0) == 1);
        ehi = elo+1;

        % Define function f = MS2-MS2cut to use for the numerical search
        MS2reps = 1e-2;
        MS2root = MS2cut-MS2reps;
        
        % Function and options for bisection search
        afun = @(TT) abs(interp_MS2_MD2_PD2(TT,eph1,eph2,HBR,Lclip) - MS2root);        
        MS2rtol = MS2reps/10;
        rbeNpts = 7;
        extrema_types = 1; % Refine minima
        endpoints = false; % Refine at endpoints
        rbeverbose = false;
        rbecheck = true;

        % Function and options for fzero search
        dfun = @(TT) interp_MS2_MD2_PD2(TT,eph1,eph2,HBR,Lclip) - MS2cut;
        fzero_options = optimset('fzero'); tolX = 1e-10;
        % fzero_options = optimset(fzero_options,'MaxIter',100);
        fzero_options = optimset(fzero_options,'Display','notify');

        % Refine the MS2 <= MS2cut segment begin times
        Ni = numel(ilo); Ti = NaN(1,Ni);
        for i=1:Ni
            % Define the root-bounding times
            Tlo = out.ref.T(ilo(i));
            Thi = out.ref.T(ihi(i));
            % Time tolerance
            Ttol = max(tolX,1e-4*(Thi-Tlo));
            % Refine the bounded root
            [TTTi,~,~,~,rbeconv,~,~,~,~,~] = ...
                refine_bounded_extrema(afun,Tlo,Thi, ...
                    rbeNpts,100,extrema_types,Ttol,MS2rtol, ...
                    endpoints,rbeverbose,rbecheck);
            if ~rbeconv || (numel(TTTi) ~= 1)
                fzero_options = optimset(fzero_options,'TolX',Ttol);
                Ti(i) = fzero(dfun,[Tlo Thi],fzero_options);
                % [Ti(i),FFi,EEi,OOi] = fzero(dfun,[Tlo Thi],fzero_options);
            else
                Ti(i) = TTTi;
            end
        end

        % Refine the MS2 <= MS2cut segment end times
        Ne = numel(elo); Te = NaN(1,Ne);
        for e=1:Ne
            % Define the root-bounding times
            Tlo = out.ref.T(elo(e));
            Thi = out.ref.T(ehi(e));
            % Time tolerance
            Ttol = max(tolX,1e-4*(Thi-Tlo));
            [TTTe,~,~,~,rbeconv,~,~,~,~,~] = ...
                refine_bounded_extrema(afun,Tlo,Thi, ...
                    rbeNpts,100,extrema_types,max(tolX,1e-4*(Thi-Tlo)),MS2rtol, ...
                    endpoints,rbeverbose,rbecheck);
            if ~rbeconv || (numel(TTTe) ~= 1)
                fzero_options = optimset(fzero_options,'TolX',Ttol);
                Te(e) = fzero(dfun,[Tlo Thi],fzero_options);
                % [Te(e),FFe,EEe,OOe] = fzero(dfun,[Tlo Thi],fzero_options);
            else
                Te(e) = TTTe;
            end
        end

        % Add the MS2 = MS2cut segment end points to the output Maha. dist.
        % reference curves
        MS2i = MS2fun(Ti); MS2e = MS2fun(Te);
        [out.ref.T,srt] = sort([out.ref.T Ti Te]);
        out.ref.MS2 = [out.ref.MS2 MS2i MS2e];
        out.ref.MS2 = out.ref.MS2(srt);
        out.ref.MD2 = [out.ref.MD2 MD2fun(Ti) MD2fun(Te)];
        out.ref.MD2 = out.ref.MD2(srt);
        if params.RefineRelDistMinima
            out.ref.RD2 = [out.ref.RD2 PD2fun(Ti) PD2fun(Te)];
            out.ref.RD2 = out.ref.RD2(srt);
        end

        % Ensure uniqueness of the time points in the reference curves
        [out.ref.T,srt] = unique(out.ref.T);
        out.ref.MS2 = out.ref.MS2(srt);
        out.ref.MD2 = out.ref.MD2(srt);
        if params.RefineRelDistMinima
            out.ref.RD2 = out.ref.RD2(srt);
        end

        % Use the reference curve times as the initial estimate of the
        % set of times of the calculated Ncdot curve
        out.T = out.ref.T;

        if debug_plotting > 0
            hold on;
            plot(Ti,MS2i,'>g');
            plot(Te,MS2e,'<y');
            hold off;
        end

    end

    % Current number of times in the Ncdot curve
    out.NT = numel(out.T);

    if debug_plotting > 0
        hold on;
        plot(out.ref.T,out.ref.MD2,':c');
        hold off;
        ymx = min([max(out.ref.MD2) (params.MScut+1)^2]);
        ymn = min(out.ref.MS2);
        ylim([ymn ymx]);
        drawnow;
    end

    % Supplement the ephemeris tables with equinoctial element states,
    % covariances and Jacobians
    if verbose
        disp(' ');
        disp('Calcuating equinoctial ephemeris quantities');
    end
    out.eph1 = calc_eq_eph(eph1,out.T(1),out.T(end));
    out.eph2 = calc_eq_eph(eph2,out.T(1),out.T(end));
    
    % Calculate Ncdot for all of the initial times
    if verbose
        disp(' ');
        disp(['Calculating ' num2str(numel(out.T)) ' initial Ncdot values']);
    end
    [Ncdot,Ucdot,Vcdot,Scdot,MDeff,CNC] = ...
        calc_Ncdot(out.T,out.eph1,out.eph2,out.ref,HBR,params); %#ok<ASGLU>
    Ncdmx = max(CNC.Rc,[],1);
    Ncdmn = min(CNC.Rc,[],1);
    
    % Curves corresponding to time points = out.T
    out.Ncdot = Ncdot;
    out.Ucdot = Ucdot;
    out.Vcdot = Vcdot;
    out.Ncdmx = Ncdmx;
    out.Ncdmn = Ncdmn;
    out.Scdot = Scdot;
    out.MDeff = MDeff;
    
    if params.SaveInitialSearch
        if verbose
            disp(' ');
            disp('Creating save mat file');
        end
        save(SaveFile,'out');
    end

end

% Calculate the Ncdot curve over the combined ephemeris times
% This is an iterative process of refining the time spacing of the Ncdot
% curve points, in order to achieve an accurate estimate of Nc, which is
% the numerical integral over the entire Ncdot curve
if verbose
    disp(' ');
    disp('Refining the Ncdot curve');
end

% Initialize refinement parameters for adaptive trapezoidal integration
refining = true; Nrefine = 1; Nnew = 0; peak_times_refined = false;
 
% Perform iterative refinement using adaptive trapezoidal integration
while refining
    
    if verbose
        disp(' ');
        disp([' Ncdot curve refinement = ' num2str(Nrefine)]);
    end
    
    if Nrefine == 1
        % Get previously calculated values for the initial grid
        Ncdot = out.Ncdot;
        Ucdot = out.Ucdot;
        Vcdot = out.Vcdot;
        Ncdmx = out.Ncdmx;
        Ncdmn = out.Ncdmn;
        Scdot = out.Scdot;
        MDeff = out.MDeff;
        NBeforeRefinement = numel(Ncdot);
    else
        % Calculate Ncdot for all new times generated during the previous
        % refinement iteration
        if ~isempty(Tnew)
            if verbose
                disp(['  Calculating ' num2str(numel(Tnew)) ' Ncdot values']);
            end
            [Ncnew,Ucnew,Vcnew,Scnew,MDnew,CNC] = ...
                calc_Ncdot(Tnew,out.eph1,out.eph2,out.ref,HBR,params); %#ok<ASGLU>
            [out.T,srt] = sort(cat(2,out.T,Tnew)); out.NT = numel(out.T);
            Ncdot = cat(2,Ncdot,Ncnew); Ncdot = Ncdot(srt);
            Ucdot = cat(2,Ucdot,Ucnew); Ucdot = Ucdot(srt);
            Vcdot = cat(2,Vcdot,Vcnew); Vcdot = Vcdot(srt);
            Scdot = cat(2,Scdot,Scnew); Scdot = Scdot(srt);
            MDeff = cat(2,MDeff,MDnew); MDeff = MDeff(srt);
            Ncdmx = cat(2,Ncdmx,max(CNC.Rc,[],1)); Ncdmx = Ncdmx(srt);
            Ncdmn = cat(2,Ncdmn,min(CNC.Rc,[],1)); Ncdmn = Ncdmn(srt);
            Nnew = Nnew + numel(Tnew);
        end
    end
    
    % Probability cutoffs
    if peak_times_refined
        probcut = probcut_iter;
    else
        probcut = probcut_last;
    end
    Nprobcut = numel(probcut);
    
    % Find all extrema in the current refinement of the Ncdot array
    [~,imxma,~,imnma] = extrema(Ncdot,true,false);
    
    % Define initial segments bracketing each Ncdot maximum
    Nseg = numel(imxma);
    % If there are no segments, break to avoid an infinite loop
	if Nseg == 0
        break;
    end
	cseg = imxma;
    iseg = NaN(size(cseg)); eseg = iseg; Prseg = iseg; Prcut = iseg;
    
    if verbose
        disp(['  Defining initial segments using ' num2str(Nseg) ' maxima']);
    end
    
    if debug_plotting == 1
        figure(2); clf;
    end

    % Define the segments evaluate statistical significance for each
    for n=1:Nseg

        % Find segment begin index
        if cseg(n) == 1
            iseg(n) = 1;
        else
            % Find preceding Ncdot minimum
            idx = find(imnma < cseg(n),1,'last');
            iseg(n) = imnma(idx);
            if (Ncdot(iseg(n)) == 0)
                % Replace Ncdot = 0 minima with MD2 maxima    
                idx = find(out.ref.TMD2max < out.T(cseg(n)),1,'last');
                TMD2mx = out.ref.TMD2max(idx);
                [~,iseg(n)] = min(abs(TMD2mx-out.T));
            end
        end
        
        % Find segment end index
        if cseg(n) == out.NT
            eseg(n) = out.NT;
        else
            % Find preceding Ncdot minimum
            idx = find(imnma > cseg(n),1,'first');
            eseg(n) = imnma(idx);
            if (Ncdot(eseg(n)) == 0)
                % Replace Ncdot = 0 minima with MD2 maxima
                idx = find(out.ref.TMD2max > out.T(cseg(n)),1,'first');
                TMD2mx = out.ref.TMD2max(idx);
                [~,eseg(n)] = min(abs(TMD2mx-out.T));
            end
        end
        
        % Determine statistical relevance of this Ncdot maximum
        cc = cseg(n);
        dd = [];
        if cseg(n) ~= iseg(n)
            dd = cat(2,dd,iseg(n));
        end
        if cseg(n) ~= eseg(n)
            dd = cat(2,dd,eseg(n));
        end
        
        % Use t-test for four samples
        Ndd = numel(dd); pp = NaN(size(dd));
        for ndd=1:Ndd
            pp(ndd) = ttest_pval(Ncdot(cc)     ,Ucdot(cc)     ,4, ...
                                 Ncdot(dd(ndd)),Ucdot(dd(ndd)),4);
        end
        [Prseg(n),mmm] = min(1-pp);
        
        % Cutoff probability for statistical relevance
        icut = min(abs(dd(mmm)-cc),Nprobcut);
        Prcut(n) = probcut(icut);
       
        if debug_plotting == 3 && Prseg(n) <= Prcut(n)
            
            % Find part of segment with positive Ncdot (padded)
            ii = iseg(n):cseg(n);
            nn = Ncdot(ii) == 0;
            if any(nn)
                ipos_n = max(ii(nn));
            else
                ipos_n = iseg(n);
            end
            ii = cseg(n):eseg(n);
            nn = Ncdot(ii) == 0;
            if any(nn)
                epos_n = min(ii(nn));
            else
                epos_n = eseg(n);
            end
            subplot(2,2,1);
            plot(out.T,Ncdot,'-k');
            hold on;
            plot(out.T,Ncdot+1.96*Ucdot,':k');
            plot(out.T,Ncdot-1.96*Ucdot,':k');
            plot(out.T,Scdot,'--g');
            j0 = cseg(n);
            j1 = iseg(n);
            j2 = eseg(n);
            plot(out.T(j0),Ncdot(j0),'^b','MarkerFaceColor','b','MarkerEdgeColor','b');
            plot(out.T([j1 j2]),Ncdot([j1 j2]),'vr','MarkerFaceColor','r','MarkerEdgeColor','r');
            idx1 = ipos_n; idx2 = epos_n;
            plot(out.T([idx1 idx2]),Ncdot([idx1 idx2]),'v', ...
                'MarkerFaceColor','y','MarkerEdgeColor','none');
            hold off;
            xrng = plot_range(out.T([j1 j2]),0.1);
            xlim(xrng);
            subplot(2,2,3);
            plot(out.ref.T,out.ref.MD2,'--k');
            hold on;
            plot(out.ref.TMD2max,out.ref.MD2max,'^b','MarkerFaceColor','b','MarkerEdgeColor','b');
            plot(out.T,MDeff,':m');
            hold off;
            xlim(xrng);
            % ndx = (out.T(j1) <= out.ref.TMD2max) & (out.ref.TMD2max <= out.T(j2));
            % yrng = plot_range(out.ref.MD2max(ndx));
            % ylim(yrng);

            subplot(2,2,[2 4]);
            plot(out.T,Ncdot,'+-k');
            hold on;
            plot(out.T,Ncdot+1.96*Ucdot,':k');
            plot(out.T,Ncdot-1.96*Ucdot,':k');
            % plot(out.T,Scdot,'--g');
            j0 = cseg(n);
            j1 = iseg(n);
            j2 = eseg(n);
            plot(out.T(j0),Ncdot(j0),'^b','MarkerFaceColor','b','MarkerEdgeColor','b');
            plot(out.T([j1 j2]),Ncdot([j1 j2]),'vr','MarkerFaceColor','r','MarkerEdgeColor','r');
            idx1 = ipos_n; idx2 = epos_n;
            plot(out.T([idx1 idx2]),Ncdot([idx1 idx2]),'v', ...
                'MarkerFaceColor','y','MarkerEdgeColor','none');
            hold off;
            xrng = plot_range(out.T([idx1 idx2]),0.1);
            xlim(xrng);
            titl = {[num2str(iseg(n)) ' ' num2str(cseg(n)) ...
                   ' ' num2str(eseg(n)) ' ' num2str(Prseg(n)) ...
                   ' ' num2str(Prcut(n))],' '};
            title(titl);
            drawnow;
            
            keyboard;
            
        end

    end
    
    % Determine the need to eliminate segments that don't have
    % statistically significant peak Ncdot values
    ndxweak = Prseg <= Prcut; % signficant probability less than cutoff
    eliminating = any(ndxweak); Nndxweak = sum(ndxweak);
    if eliminating
        if verbose
            disp(['  Statistically weak Ncdot maxima = ' num2str(Nndxweak)]);
        end
        if params.RetainOneSegment && Nseg == 1
            eliminating = false;
            if verbose
                disp('  Retaining least statistically weak Ncdot maximum');
            end
        end
    end
    
    % Eliminate all segments that don't have statistically significant
    % peak Ncdot values, i.e., those with weak maxima. Eliminate the
    % weakest first, and iteratively update the evaluated weaknesses
    while eliminating
        
        % Sort the weak maxima into increasing order
        ndx = find(ndxweak);
        [~,srt] = sort(Ncdot(cseg(ndx)));
        ndx = ndx(srt);
        
        % Select the lowest to eliminate
        n = ndx(1);
        
        % Allocate new segment center, initial and ending indices after
        % elimination
        cnew = cseg; inew = iseg; enew = eseg;
        
        % Mark current segment for deletion using an NaN value
        cnew(n) = NaN;

        % Determine if this is a standalone maximum, or if it abuts another        
        abutlft = any(eseg == iseg(n));
        abutrgt = any(iseg == eseg(n));
        
        % Process abutting segments
        if abutlft || abutrgt
            
            % Decide which of the two bracketing segments to eliminate
            if iseg(n) == 1
                elim = eseg(n);
            elseif eseg(n) == out.NT
                elim = iseg(n);
            else
                % If abutted on both sides, eliminate the bracketing
                % minimum point that has the smallest statistical
                % significance
                if abutlft && abutrgt
                    % Statistical z-value for init point (min. on left)
                    qq = Ncdot(cseg(n))-Ncdot(iseg(n));
                    dqq = sqrt(Ucdot(cseg(n))^2 + Ucdot(iseg(n))^2);
                    izz2 = (qq/dqq)^2;
                    % Statistical z-value for end point (min. on right)
                    qq = Ncdot(cseg(n))-Ncdot(eseg(n));
                    dqq = sqrt(Ucdot(cseg(n))^2 + Ucdot(eseg(n))^2);
                    ezz2 = (qq/dqq)^2;
                    % Eliminate min point w/ smaller statistical z value
                    if izz2 < ezz2
                        elim = iseg(n);
                    else
                        elim = eseg(n);
                    end
                elseif abutlft
                    elim = iseg(n);
                else % if abutrght
                    elim = eseg(n);
                end
            end
            
            % Eliminate the weakest max and associated min, and extend the
            % adjacent segment if abutted
            if elim == iseg(n)
                m = n-1;
                if enew(m) == iseg(n)
                    enew(m) = eseg(n);
                end
            else
                m = n+1;
                if inew(m) == eseg(n)
                    inew(m) = iseg(n);
                end
            end
            
        end

        % Debug plotting to see elimination process
        if debug_plotting == 2 && Ncdot(cseg(n)) > 5e-8

            subplot(1,2,1);
            plot(out.T,Ncdot,'.-k');
            hold on;
            plot(out.T,Ncdot+1.96*Ucdot,':k');
            plot(out.T,Ncdot-1.96*Ucdot,':k');
            % plot(out.T,Scdot,'--g');
            jmin = 1; jmax = Nseg;
            % j0 = unique([cseg(max(jmin,n-1)) cseg(n) cseg(min(jmax,n+1))]);
            j1 = [iseg(max(jmin,n-2)) iseg(max(jmin,n-1)) iseg(n) iseg(min(jmax,n+1)) iseg(min(jmax,n+2))];
            j2 = [eseg(max(jmin,n-2)) eseg(max(jmin,n-1)) eseg(n) eseg(min(jmax,n+1)) eseg(min(jmax,n+2))];
            j12 = unique([j1 j2]);
            % plot(out.T(j12),Ncdot(j12),'v','MarkerFaceColor','m','MarkerEdgeColor','m');%
            % plot(out.T(j0),Ncdot(j0),'^','MarkerFaceColor','c','MarkerEdgeColor','c');
            ieseg = unique([iseg eseg]);
            plot(out.T(ieseg),Ncdot(ieseg),'v','MarkerFaceColor','m','MarkerEdgeColor','m');
            plot(out.T(cseg),Ncdot(cseg),'^','MarkerFaceColor','c','MarkerEdgeColor','c');
            hold off;
            xrng = plot_range(out.T(j12),0.1);
            xlim(xrng);
            titl = {[num2str(iseg(n)) ' ' num2str(cseg(n)) ...
                   ' ' num2str(eseg(n)) ' ' num2str(Prseg(n)) ...
                   ' ' num2str(Prcut(n))],' '};
            title(titl);

            subplot(1,2,2);
            plot(out.T,Ncdot,'.-k');
            hold on;
            plot(out.T,Ncdot+1.96*Ucdot,':k');
            plot(out.T,Ncdot-1.96*Ucdot,':k');
            % plot(out.T,Scdot,'--g');
            if abutlft || abutrgt
                j0 = cnew(m);
                j1 = max(1,inew(m)-1);
                j2 = min(numel(out.T),enew(m)+1);
                plot(out.T([j1 j2]),Ncdot([j1 j2]),'v','MarkerFaceColor','r','MarkerEdgeColor','r');
                plot(out.T(j0),Ncdot(j0),'^','MarkerFaceColor','b','MarkerEdgeColor','b');
                titl = {[num2str(inew(m)) ' ' num2str(cnew(m)) ...
                       ' ' num2str(enew(m))],' '};
                title(titl);
            end
            hold off;
            xrng = plot_range(out.T(unique([j12 j1 j2])),0.1);
            xlim(xrng);

            drawnow;

            keyboard;

        end
        
        % Define the new segments
        ndx = ~isnan(cnew);
        cseg = cnew(ndx); 
        iseg = inew(ndx);
        eseg = enew(ndx);
        Nseg = numel(cseg);

        % Update the statistical relevance of uneliminated Ncdot maxima
        Prseg = NaN(size(cseg)); Prcut = Prseg;        
        for n=1:Nseg
            % Determine statistical relevance of maximum
            cc = cseg(n);
            dd = [];
            if cseg(n) ~= iseg(n)
                dd = cat(2,dd,iseg(n));
            end
            if cseg(n) ~= eseg(n)
                dd = cat(2,dd,eseg(n));
            end
            % % Estimate updated z-value, probability and prob. cutoff
            % qq = Ncdot(cc)-Ncdot(dd);
            % dqq = sqrt(Ucdot(dd).^2 + Ucdot(cc).^2);
            % zz = qq./dqq; zz2 = zz.^2;
            % [Prseg(n),mmm] = min(chi2cdf(zz2,1));
            % Use t-test for four samples
            Ndd = numel(dd); pp = NaN(size(dd));
            for ndd=1:Ndd
                pp(ndd) = ttest_pval(Ncdot(cc)     ,Ucdot(cc)     ,4, ...
                                     Ncdot(dd(ndd)),Ucdot(dd(ndd)),4);
            end
            [Prseg(n),mmm] = min(1-pp);
            % Cutoff probability
            icut = min(abs(dd(mmm)-cc),Nprobcut);
            Prcut(n) = probcut(icut);
        end
        
        % Continue eliminating if there are still any weak Ncdot maxima
        ndxweak = Prseg <= Prcut;
        eliminating = any(ndxweak);
        if params.RetainOneSegment && Nseg == 1
            eliminating = false;
            if verbose
                disp('  Retaining the least statistically weak Ncdot maximum');
            end
        end
        
    end
    
    % Process the refined segments (if any)
    if Nseg > 0

        % Perform segment-by-segment trapezoidal integration refinement of
        % the Ncdot curve time steps
        if verbose
            disp(['  Refining time steps for retained time segments = ' ...
                num2str(Nseg)]);
        end

        % if debug_plotting > 1
        %     figure(2); clf;
        % end

        ndxrefine = [];

        Ncseg = NaN(size(iseg)); Ucseg = Ncseg;
        Ncsmx = Ncseg; Ncsmn = Ncseg;
        ipos  = iseg; epos  = eseg;

        % Segment-by-segement trapezoidal integration
        for n=1:Nseg

            % Find part of segment with positive Ncdot (padded)
            ii = iseg(n):cseg(n);
            nn = Ncdot(ii) == 0;
            if any(nn)
                ipos(n) = max(ii(nn));
            end
            ii = cseg(n):eseg(n);
            nn = Ncdot(ii) == 0;
            if any(nn)
                epos(n) = min(ii(nn));
            end

            % Calculate the segment Nc
            idx = ipos(n):epos(n); Nidx = numel(idx); Midx = Nidx-1;
            Tseg = out.T(idx);
            Tdif = diff(Tseg);
            dT = [Tdif(1) Tdif(1:Midx-1)+Tdif(2:Midx) Tdif(Midx)]/2;
            dTs = dT*86400;
            NNseg = dTs.*Ncdot(idx);
            Ncseg(n) = sum(NNseg);
            Ncsmx(n) = sum(dTs.*Ncdmx(idx));
            Ncsmn(n) = sum(dTs.*Ncdmn(idx));
            
            % Calculate the segment Nc uncertainty
            UUseg = dTs.*Ucdot(idx);
            Ucseg(n) = sqrt(sum(UUseg.*UUseg));
            if Ncseg(n) > 0
                NcdotAvg = sum(NNseg.*Ncdot(idx))/Ncseg(n);
                UcdotAvg = sum(NNseg.*Ucdot(idx))/Ncseg(n);
                Ucseg(n) = max([Ucseg(n), ...
                    Ncseg(n)*UcdotAvg/NcdotAvg]);
            end

            % Find trapezoidal time steps that need bisection refinement
            if Nrefine == 1
                % First refinement:
                % Ensure steps have sufficiently small trapezoidal integ
                % components, and sufficiently small time steps
                ndx = (dT > Tcsegtol*(Tseg(Nidx)-Tseg(1))) | ...
                      (NNseg > Ncsegtol*Ncseg(n));
                if isempty(ndx)
                    % Always refine at least one initially
                    ndx = NNseg == max(NNseg);
                end
            else
                % Ensure steps have sufficiently small trapezoidal integ
                % components for both this segment and the total over segments
                if (Ncseg(n) > params.NcFracTiny*Nctot)
                    ndx = (dT > Tcsegtol*(Tseg(Nidx)-Tseg(1))) | ...
                          (NNseg > Ncsegtol*Ncseg(n)) | ...
                          (NNseg > Nctottol*Nctot);
                elseif (Ncseg(n) > params.NcFracNegligible*Nctot)
                    ndx = (dT > Tctnytol*(Tseg(Nidx)-Tseg(1))) | ...
                          (NNseg > Nctnytol*Ncseg(n));
                else
                    ndx = [];
                end
            end

            % Add bisection refinements for this segment to total list of
            % refinements
            if any(ndx)
                ndxrefine = cat(2,ndxrefine,idx(ndx));
            end

            if debug_plotting == 1
                
                figure(2); clf;

                subplot(2,2,1);
                plot(out.T,Ncdot,'-k');
                hold on;
                plot(out.T,Ncdot+1.96*Ucdot,':k');
                plot(out.T,Ncdot-1.96*Ucdot,':k');
                plot(out.T,Scdot,'--g');
                j0 = cseg(n);
                j1 = iseg(n);
                j2 = eseg(n);
                plot(out.T(j0),Ncdot(j0),'^b','MarkerFaceColor','b','MarkerEdgeColor','b');
                plot(out.T([j1 j2]),Ncdot([j1 j2]),'vr','MarkerFaceColor','r','MarkerEdgeColor','r');
                idx1 = ipos(n); idx2 = epos(n);
                plot(out.T([idx1 idx2]),Ncdot([idx1 idx2]),'v', ...
                    'MarkerFaceColor','y','MarkerEdgeColor','none');
                hold off;
                xrng = plot_range(out.T([j1 j2]),0.1);
                xlim(xrng);
                subplot(2,2,3);
                semilogy(NaN,NaN); hold on;
                plot(out.T,MDeff,'.-m','LineWidth',1);
                plot(out.ref.T,out.ref.MD2,'-k','LineWidth',0.5);
                plot(out.ref.TMD2max,out.ref.MD2max,'^b','MarkerFaceColor','b','MarkerEdgeColor','b');
                hold off;
                xlim(xrng);

                subplot(2,2,[2 4]);
                plot(out.T,Ncdot,'+-k');
                hold on;
                plot(out.T,Ncdot+1.96*Ucdot,':k');
                plot(out.T,Ncdot-1.96*Ucdot,':k');
                % plot(out.T,Scdot,'--g');
                j0 = cseg(n);
                j1 = iseg(n);
                j2 = eseg(n);
                plot(out.T(j0),Ncdot(j0),'^b','MarkerFaceColor','b','MarkerEdgeColor','b');
                plot(out.T([j1 j2]),Ncdot([j1 j2]),'vr','MarkerFaceColor','r','MarkerEdgeColor','r');
                idx1 = ipos(n); idx2 = epos(n);
                plot(out.T([idx1 idx2]),Ncdot([idx1 idx2]),'v', ...
                    'MarkerFaceColor','y','MarkerEdgeColor','none');
                hold off;
                xrng = plot_range(out.T([idx1 idx2]),0.1);
                xlim(xrng);
                titl = [num2str(ipos(n)) ' ' num2str(cseg(n)) ...
                       ' ' num2str(epos(n))];
                title(titl);
                drawnow;

                keyboard;

            end

        end

        % Sum over all segments
        Nctot = sum(Ncseg);
        Nctmx = sum(Ncsmx);
        Nctmn = sum(Ncsmn);
        
        if verbose
            Uctot = sqrt(sum(Ucseg.^2));
            [~,estr] = smart_error_format(Nctot,1.96*Uctot);
            [xmdstr,xlostr,xhistr] = smart_error_range( ...
                Nctot,Nctot-1.96*Uctot,Nctot+1.96*Uctot);
            disp(['  Nctot = ' xmdstr ' +/- ' estr ...
                  ' (95% ' xlostr ' to ' xhistr ')']);
            [~,xlostr,xhistr] = smart_error_range( ...
                Nctot,Nctmn,Nctmx);
            disp(['  Nctot Unc. Range = ' xlostr ' to ' xhistr]);
            upct = 100*max((Nctmx-Nctot)/Nctot,(Nctot-Nctmn)/Nctot);
            disp(['  Nctot Unc. Ratio = ' smart_exp_format(upct,3) '%']);
        end

        % Perform integ time step refinements, to ensure each segment has
        % adequate integ spacing, and no time steps are too large w.r.t.
        % the total of the integration.  As a last step, refine the peak
        % times using parabolic fitting
        if isempty(ndxrefine)
            if peak_times_refined
                % Refining complete
                refining = false;
                if verbose
                    disp(['  ' num2str(Nnew) ' integration times added after ' ...
                        num2str(Nrefine) ' refinements']);
                end    
            else
                % One more process step to refine the peak times
                peak_times_refined = true;
                Tpeak = NaN(size(cseg));
                d = 2; % Points on either side of peak used for parabolic fit
                Nseg = numel(cseg);
                for n=1:Nseg
                    c = cseg(n);
                    if ~( (c == 1) || (c == out.NT) )
                        m = max(1,c-d); p = min(out.NT,c+d); i = m:p;
                        dT = out.T(p)-out.T(m);
                        dN = max(Ncdot(i))-min(Ncdot(i));
                        x = (out.T(i)-out.T(c))/dT;
                        y = (Ncdot(i)-Ncdot(c))/dN;
                        [P,~] = polyfit(x,y,2);
                        x0 = -P(2)/2/P(1);
                        y0 = P(1)*x0^2+P(2)*x0+P(3);
                        % Only accept the the parabolic fit time if that y
                        % value is higher than the y value of the fitted
                        % points
                        if y0 >= max(y)
                            Tpeak(n) = out.T(c)+dT*x0;
                        end
                    end
                end
                % Incorporate the refined peak times
                ndx = ~isnan(Tpeak) & ~ismember(Tpeak,out.T);
                Tnew = unique(Tpeak(ndx));
                if verbose
                    disp(['  ' num2str(numel(Tnew)) ' peak times fit after ' ...
                        num2str(Nrefine) ' refinements']);
                end   
            end
        else
            % Define the set of time steps to refine on the next iteration
            Tnew = [];
            Nndxrefine = numel(ndxrefine);
            for rr=1:Nndxrefine
                r = ndxrefine(rr);
                if r > 1
                    % Bisect earlier time step
                    Tnew = cat(2,Tnew,0.5*(out.T(r-1)+out.T(r)));
                end
                if r < out.NT
                    % Bisect later time step
                    Tnew = cat(2,Tnew,0.5*(out.T(r)+out.T(r+1)));
                end
            end
            Tnew = unique(Tnew);
            % Increment refinement counter
            Nrefine = Nrefine+1;
            if verbose
                disp(['  Number of new integration times = ' num2str(numel(Tnew))]);
            end
        end

        if debug_plotting > 0

            figure(3); clf;

            subplot(2,1,1);
            plot(out.T,Ncdot,'*-k');
            hold on;
            plot(out.T,Ncdot-1.96*Ucdot,':k');
            plot(out.T,Ncdot+1.96*Ucdot,':k');
            plot(out.T(cseg),Ncdot(cseg),'^b');
            plot(out.T(iseg),Ncdot(iseg),'vr');
            plot(out.T(eseg),Ncdot(eseg),'vr');
            plot(out.T(ipos),Ncdot(ipos),'v', ...
                'MarkerFaceColor','y','MarkerEdgeColor','none');
            plot(out.T(epos),Ncdot(epos),'v', ...
                'MarkerFaceColor','y','MarkerEdgeColor','none');
            hold off;
            xrng = plot_range(out.T,0.1);
            xlim(xrng);
            subplot(2,1,2);
            plot(out.ref.T,out.ref.MD2,'--k');
            hold on;
            plot(out.ref.TMD2max,out.ref.MD2max,'^b');
            plot(out.T,MDeff,':m');
            hold off;
            xlim(xrng);

            % keyboard;

        end
        
    end
    
end

if verbose
    NAfterRefinement = numel(Ncdot); 
    disp(' ');
    disp(['Before refinement: ' num2str(NBeforeRefinement) ' Ncdot values']);
    disp([' After refinement: ' num2str(NAfterRefinement ) ' Ncdot values']);
end

% Assemble output quantities

% Curves corresponding to time points = out.T
out.Ncdot = Ncdot;
out.Ucdot = Ucdot;
out.Vcdot = Vcdot;
out.Ncdmx = Ncdmx;
out.Ncdmn = Ncdmn;
out.Scdot = Scdot;
out.MDeff = MDeff;

out.Nseg  = Nseg;
if Nseg > 0
    out.cseg  = cseg;  
    out.iseg  = iseg;
    out.eseg  = eseg;
    out.ipos  = ipos; 
    out.epos  = epos;
    out.Ncseg = Ncseg;
    out.Ucseg = Ucseg;
    out.Ncsmx = Ncsmx;
    out.Ncsmn = Ncsmn;
end

% Cumulative Nc curve and interpolation uncertainty curve
Tdif = diff(out.T)*86400;
dT = [Tdif(1) Tdif(1:end-1)+Tdif(2:end) Tdif(end)]/2;
out.Nccum = cumsum(dT.*Ncdot);
out.Uccum = sqrt(cumsum((dT.*Ucdot).^2));
out.Nccmx = cumsum(dT.*Ncdmx);
out.Nccmn = cumsum(dT.*Ncdmn);

% Final Nc value and interpolation uncertainty
if Nseg > 0
    out.Nc = sum(out.Ncseg);
    out.Uc = sqrt(sum(out.Ucseg.^2));
    out.Ncmx = sum(out.Ncsmx);
    out.Ncmn = sum(out.Ncsmn);
else
    out.Nc = out.Nccum(end); 
    out.Uc = out.Uccum(end);
    out.Ncmx = out.Nccmx(end);
    out.Ncmn = out.Nccmn(end);
end
Nc = out.Nc;

% Calculate the estimated cumulative Pc, and the min/max Pc bounds

% Initialize the min/max Pcdot curves
PcdotMin = out.Ncdot;
PcdotMax = out.Ncdot;

% Time differences in seconds
Tdif = diff(out.T)*86400;
dT = [Tdif(1) Tdif(1:end-1)+Tdif(2:end) Tdif(end)]/2;

% Initialize cumulative upper and lower limit Pc values
if out.Nseg > 0    
    PcumMaxLast = out.Ncseg(1); 
    PcumMinLast = out.Ncseg(1);
else
    % If there are no distinct time segments with separate collision rate 
    % peaks, use the total collision probability over the entire interval, 
    % which is represented in out.Nc
    PcumMaxLast = out.Nc;
    PcumMinLast = out.Nc;
end

if out.Nseg > 1

    for n=2:out.Nseg

        % Min and max cumulative probs
        PcumMaxNext = min( 1 , -out.Ncseg(n)*PcumMaxLast + out.Ncseg(n) + PcumMaxLast );            
        PcumMinNext = max( PcumMinLast, out.Ncseg(n) );

        % Ncdot > 0 begin and end time indices for segment
        i = out.ipos(n); e = out.epos(n);

        % Adjust segment Pcdot values to sum to Pcum max limit
        ratio = (PcumMaxNext-PcumMaxLast)/out.Ncseg(n);
        PcdotMax(i:e) = out.Ncdot(i:e)*ratio;

        % Adjust segment Pcdot values to sum to PcumMin limit
        ratio = (PcumMinNext-PcumMinLast)/out.Ncseg(n);
        PcdotMin(i:e) = out.Ncdot(i:e)*ratio;

        % Update cumulative probs
        PcumMaxLast = PcumMaxNext;
        PcumMinLast = PcumMinNext;

    end

end

% Last Pcum values represent min/max Pc estimates
out.PcMax = PcumMaxLast; PcMax = out.PcMax;
out.PcMin = PcumMinLast; PcMin = out.PcMin;

% Calculate min/max cumulative Pc curves
out.PcumMax = cumsum(dT.*PcdotMax);
out.PcumMin = cumsum(dT.*PcdotMin);

% Create save file
if params.MakeOutFile
    OutFile = [fullfile(params.outputPath,params.outputRoot) '.mat'];
    save(OutFile,'out');
end

return;
end

% =========================================================================

function [Ncdot,Ucdot,Vcdot,Scdot,MDeff,out] = calc_Ncdot(T,eph1,eph2,ref,HBR,params)

% Calculate Ncdot values for the times in array T, along with the
% subcomponent rates Rc and associated weights Wc

% Extract bracketing flag
bracketing = params.bracketing;

% Extract UseEqEls flag
UseEqEls = params.UseEqEls;

% Constants
twopi = 2*pi; twopicubed = twopi^3;
H = HBR; H2 = H^2;

% Eigenvalue clipping parameter (km^2)
Lclip = (params.Fclip*HBR)^2;

% Number of time points
NT = numel(T);

% Allocate main output arrays
Ncdot = NaN(size(T)); 
Ucdot = Ncdot; 
Vcdot = Ncdot;
Scdot = Ncdot;
MDeff = Ncdot;

% Set up the rate arrays and weights
if ~bracketing
    % Nearest neighbor only has one weight
    NW = 1; % rNW = 1;
else
    % Bracketing eph points result in four combininations and weights
    % NW = 4; % rNW = NW/(NW-1);
    NW = 5; % rNW = NW/(NW-1);
end
out.Rc = NaN(NW,NT);
out.Wc = NaN(NW,NT);
    
% Interpolate MS2 to determine whichs rates can be inferred as zero
MS2 = interp1(ref.T,ref.MS2,T);
ndx0 = MS2 > params.MSZeroNcdot^2;
out.Rc(:,ndx0) = 0; out.Wc(:,ndx0) = 1/NW;
Ncdot(ndx0) = 0; Ucdot(ndx0) = 0; Vcdot(ndx0) = 0; MDeff(ndx0) = Inf;
ndx = find(~ndx0);
Nndx = numel(ndx);

% Variables for parallel execution
tt = T(ndx);

Ncdt = NaN(Nndx,1);
Ucdt = NaN(Nndx,1);
Vcdt = NaN(Nndx,1);
Scdt = NaN(Nndx,1);
MDef = NaN(Nndx,1);

Rc = NaN(Nndx,NW);
Wc = NaN(Nndx,NW);
% Sc = NaN(Nndx,NW);
% MD = NaN(Nndx,NW);

params_use_Lebedev = params.use_Lebedev;
params_vec_Lebedev = params.vec_Lebedev;
params_wgt_Lebedev = params.wgt_Lebedev;
if params_use_Lebedev
    params_AbsTol = [];
    params_RelTol = [];
    params_MaxFunEvals = [];
else
    params_AbsTol = 0;
    params_RelTol = 1e-9;
    params_MaxFunEvals = 2000;    
end

% Initialize peak overlap determination parameters
POPEpar0.maxiter = params.POPmaxiter;

% Calculate rates for all times with MS2 <= MS2cut
parfor nnT=1:Nndx % using_parfor
    
    % Current time
    t = tt(nnT);
    
    % Calculate states, covs and Jacobians from bracketing points
    [c1,X1c,P1c,E1c,Q1c,J1c,W1c,b1,X1b,P1b,E1b,Q1b,J1b,W1b] = ...
        BracketPoints(t,eph1,bracketing,UseEqEls); %#ok<ASGLU>
    [c2,X2c,P2c,E2c,Q2c,J2c,W2c,b2,X2b,P2b,E2b,Q2b,J2b,W2b] = ...
        BracketPoints(t,eph2,bracketing,UseEqEls); %#ok<ASGLU>

    % Allocate arrays for weighted uncertainty estimates
    WW = NaN(1,NW);
    RR = NaN(1,NW);
    SS = NaN(1,NW);
    MM = NaN(1,NW);

    X1 = []; P1 = []; J1 = []; E1 = []; Q1 = []; W1 = [];
    X2 = []; P2 = []; J2 = []; E2 = []; Q2 = []; W2 = []; %#ok<NASGU>
    
    for nW=1:NW
        
        % States and Covariances
        if nW == 1
            X1 = X1c; P1 = P1c; E1 = E1c; Q1 = Q1c; J1 = J1c; W1 = W1c;
            X2 = X2c; P2 = P2c; E2 = E2c; Q2 = Q2c; J2 = J2c; W2 = W2c;
        elseif nW == 2
           %X1 = X1c; P1 = P1c; E1 = E1c; Q1 = Q1c; J1 = J1c; W1 = W1c;
            X2 = X2b; P2 = P2b; E2 = E2b; Q2 = Q2b; J2 = J2b; W2 = W2b;
        elseif nW == 3
            X1 = X1b; P1 = P1b; E1 = E1b; Q1 = Q1b; J1 = J1b; W1 = W1b;
            X2 = X2c; P2 = P2c; E2 = E2c; Q2 = Q2c; J2 = J2c; W2 = W2c;
        elseif nW == 4
           %X1 = X1b; P1 = P1b; E1 = E1b; Q1 = Q1b; J1 = J1b; W1 = W1b;
            X2 = X2b; P2 = P2b; E2 = E2b; Q2 = Q2b; J2 = J2b; W2 = W2b;
        else % if nW == 5 - Blended state and cov
            
            % No weights for the blended value
            W1 = 0; W2 = 0;
            % Blended PV states and covs
            X1 = W1c*X1c + W1b*X1b;
            P1 = W1c*P1c + W1b*P1b;
            X2 = W2c*X2c + W2b*X2b;
            P2 = W2c*P2c + W2b*P2b;
            % Calculate the equinoctial elements and Jacobians
            if UseEqEls
                [~,n,af,ag,chi,psi,lM,~] = ...
                    convert_cartesian_to_equinoctial(X1(1:3),X1(4:6));
                E1 = [n; af; ag; chi; psi; lM];
                J1 = jacobian_equinoctial_to_cartesian(E1,X1);
                I6x6 = eye(6,6);
                K1 = J1\I6x6;
                Q1 = K1 * P1 * K1';
                [~,n,af,ag,chi,psi,lM,~] = ...
                    convert_cartesian_to_equinoctial(X2(1:3),X2(4:6));
                E2 = [n; af; ag; chi; psi; lM];
                J2 = jacobian_equinoctial_to_cartesian(E2,X2);
                K2 = J2\I6x6;
                Q2 = K2 * P2 * K2';
            else
                % Use best estimate state, and calculate eph. offset quantities
                [E1,C2E] = convert_cartesian_to_ephoffset(t,X1,X1,eph1);
                M3 = [C2E.Xhat_tau  ; C2E.Yhat_tau  ; C2E.Zhat_tau];
                Z3 = zeros(3,3);
                J1 = [M3 Z3; Z3 M3]; % dE/dX
                Q1 = J1 * P1 * J1';
                [E2,C2E] = convert_cartesian_to_ephoffset(t,X2,X2,eph2);
                M3 = [C2E.Xhat_tau  ; C2E.Yhat_tau  ; C2E.Zhat_tau];
                J2 = [M3 Z3; Z3 M3]; % dE/dX
                Q2 = J2 * P2 * J2';
            end
            
        end
        
        % Weight
        WW(nW) = W1*W2;

        % Find the peak PD overlap position
        
        % Initialize POP parameters
        POPEpar = POPEpar0;

        if UseEqEls
            % Use POP function in equinoctial element mode
            POPEpar.Jb1 = J1;
            POPEpar.Eb1 = E1;
            POPEpar.Qb1 = Q1;
            POPEpar.Jb2 = J2;
            POPEpar.Eb2 = E2;
            POPEpar.Qb2 = Q2;
            [converged,~,~,~,POP] = PeakOverlapPosEph( ...
                t,X1,P1,X2,P2,[],[],HBR,POPEpar);
        else
            % Use POP function in eph. offset element mode
            [converged,~,~,~,POP] = PeakOverlapPosEph( ...
                t,X1,P1,X2,P2,eph1,eph2,HBR,POPEpar);
        end
        
        % If the POP function converged, then calculate Rc = Ncdot values

        if converged

            % Calculate the relative offset pos and vel values
            ru = POP.xu2(1:3)-POP.xu1(1:3); 
            vu = POP.xu2(4:6)-POP.xu1(4:6);

            % Calculate pos/vel covariances
            Ps1 = cov_make_symmetric(POP.Js1 * Q1 * POP.Js1');
            Ps2 = cov_make_symmetric(POP.Js2 * Q2 * POP.Js2');

            % Calculate the relative pos/vel coviariance, extract the
            % submatrices, and calculate related quantities
            Ps = Ps1+Ps2; As = Ps(1:3,1:3); Bs = Ps(4:6,1:3); Cs = Ps(4:6,4:6);
            [~,~,~,~,~,Asdet,Asinv] = CovRemEigValClip(As,Lclip);
            % Ns0 = (twopicubed*Asdet)^(-0.5);
            Ns0 = 1/sqrt(twopicubed*Asdet);
            bs = Bs*Asinv; Csp = Cs-bs*Bs';

            % Calculate the effective Mahalanobis distance squared
            MM(nW) = ru' * Asinv * ru;
            
            % Calculate the small-HBR limit coll. rate approximation
            Ns = Ns0 * exp(-0.5*MM(nW));
            SS(nW) = H2 * pi * Ns * norm(vu);

            % Calculate the Ncdot value using numerical integration, over
            % the unit sphere, valid for any HBR value

            if params_use_Lebedev

                % Use Lebedev for integration over the unit sphere
                Pint = Ncdot_integrand(params_vec_Lebedev,ru,vu,Asinv, ...
                    H,bs,Csp,-log(Ns0),false);
                szPint = size(Pint);
                if (szPint(1) ~= 1); Pint = Pint'; end
                RR(nW) = H2 * Pint * params_wgt_Lebedev;

            else

                % Define the anonymous function for the quad2d integrand
                fun = @(ph,u)Ncdot_quad2d_integrand(ph,u, ...
                    ru,vu,Asinv,H,bs,Csp,-log(Ns0),false);
                % Perform the quad2d integration
                [Pint,~] = quad2d(fun,0,twopi,-1,1, ...
                    'AbsTol',params_AbsTol,'RelTol',params_RelTol, ...
                    'MaxFunEvals',params_MaxFunEvals);
                RR(nW) = H2 * Pint;
                
            end
            
        else
            
            % If not converged, and failed because negative energy states
            % were achieved, assume zero rates
            if POP.failure ~= 0
                SS(nW) = 0;
                RR(nW) = 0;
            else
                warning('POP determinination unconverged');
                SS(nW) = 0;
                RR(nW) = 0;
            end
            
        end
   
    end
    
    if abs(1-sum(WW)) > 10*eps
        warning('WW sum bad');
    end

    Rc(nnT,:) = RR;
    Wc(nnT,:) = WW;
    % Sc(nnT,:) = SS;
    % MD(nnT,:) = MM;
    
    Ncdt(nnT) = RR(5);
    Scdt(nnT) = SS(5);
    MDef(nnT) = MM(5);
    
    if bracketing
        Vcdt(nnT) = std(RR(1:4),WW(1:4));
        Ucdt(nnT) = std(RR(1:4));
    else
        Ucdt(nnT) = 0; Vcdt(nnT) = 0;
    end

end

ndx0 = ~ndx0;
out.Rc(:,ndx0) = Rc';
out.Wc(:,ndx0) = Wc';

Ncdot(ndx0) = Ncdt';
Ucdot(ndx0) = Ucdt';
Vcdot(ndx0) = Vcdt';

Scdot(ndx0) = Scdt';
MDeff(ndx0) = MDef';

return
end

% =========================================================================

function [c,Xc,Pc,Ec,Qc,Jc,Wc,b,Xb,Pb,Eb,Qb,Jb,Wb] = BracketPoints(t,eph,bracketing,UseEqEls)

% Calculate the PV state X and Jacobian matrix dX/dE at time t using
% the closest ephemeris points and (optionally) the bracketing point

% Find the bracketing indices
[c,b] = find_bracketing_indices(t,eph.T);

% Weight by bracketing time differences
if bracketing
    Wb = (t-eph.T(c)) / (eph.T(b)-eph.T(c));
    if (Wb < 0) || (Wb > 1)
        error('Invalid bracketing weight');
    end 
    Wc = 1 - Wb;
else
    Wc = 1; Wb = 0;
end

% Calculate the state interpolants from closest and bracketing points
[ri,vi,~,rc,vc,~,rb,vb] = interpState(t,eph);

% Calculate the covariance interpolant from the closest point
rc = rc'; vc = vc'; Xc = [rc; vc];
[~,Pc,~,~,Qc,~] = TBSTMCovInterp(t,Xc, ...
    eph.T(c),eph.X(:,c),eph.P(:,:,c),  ...
    eph.T(b),eph.X(:,b),eph.P(:,:,b));

% Calculate the equinoctial elements and Jacobian for the closest point
if UseEqEls
    [~,n,af,ag,chi,psi,lM,~] = ...
        convert_cartesian_to_equinoctial(rc,vc);
    Ec = [n; af; ag; chi; psi; lM];
    Jc = jacobian_equinoctial_to_cartesian(Ec,Xc);
else
    % Use best estimate state, and calculate eph. offset quantities
    Xc = [ri vi]';
    [Ec,C2E] = convert_cartesian_to_ephoffset(t,Xc,Xc,eph);
    M3 = [C2E.Xhat_tau  ; C2E.Yhat_tau  ; C2E.Zhat_tau];
    Z3 = zeros(3,3);
    Jc = [M3 Z3; Z3 M3]; % dE/dX
    Qc = Jc * Pc * Jc';
end

% Calculate the covariance interpolant for the bracketing point
rb = rb'; vb = vb'; Xb = [rb; vb];
[~,~,Pb,~,~,Qb] = TBSTMCovInterp(t,Xb, ...
    eph.T(c),eph.X(:,c),eph.P(:,:,c),  ...
    eph.T(b),eph.X(:,b),eph.P(:,:,b));

% Calculate the equinoctial elements and Jacobian for the
% bracketing point
if UseEqEls
    [~,n,af,ag,chi,psi,lM,~] = ...
        convert_cartesian_to_equinoctial(rb,vb);
    Eb = [n; af; ag; chi; psi; lM];
    Jb = jacobian_equinoctial_to_cartesian(Eb,Xb);
else
    % Use best estimate state, and calculate eph. offset quantities
    Xb = Xc; Eb = Ec; Jb = Jc;
    Qb = Jb * Pb * Jb';    
end

return
end

% =========================================================================

function out = calc_eq_eph(eph,taua,taub)

% Limit an ephemeris to span the interval taua to taub, and supplement with
% equinoctial states, covariances, and Jacobians if necessary

% Find indices spanning (taua,taub)
ndx = find((taua <= eph.T) & (eph.T <= taub));
if isempty(ndx)
    [~,na] = min(abs(taua-eph.T));
    [~,nb] = min(abs(taub-eph.T));
    ndx = unique([na nb]);
end

% Pad by three to allow 5-point state interpolation
na = max(ndx(1)-3,1);
nb = min(ndx(end)+3,eph.N);
ndx = na:nb;

% Copy limited times, PV states and PV covariances
out.T = eph.T(ndx);
out.N = numel(out.T);
out.X = eph.X(:,ndx);
out.P = eph.P(:,:,ndx);

% Check if the equinoctial states, covariances and Jacobians exist

if isfield(eph,'E') % assumes Q, J, and K exist too
    
    % Copy the equinoctial states, covariances and Jacobians
    out.E = eph.E(:,ndx);
    out.Q = eph.Q(:,:,ndx);
    out.J = eph.J(:,:,ndx);
    out.K = eph.K(:,:,ndx);

else

    % Allocate the equinoctial states, covariances and Jacobians
    out.E = NaN(6,out.N);
    out.Q = NaN(6,6,out.N);
    out.J = NaN(6,6,out.N);
    out.K = NaN(6,6,out.N);

    % Generate the equinoctial states, covariances and Jacobians
    
    I6x6 = eye(6,6); twopi = 2*pi;

    for n=1:out.N

        % Generate the equinoctial elements from the PV vector
        [~,nmean,afmean,agmean,chimean,psimean,lMmean,~] = ...
            convert_cartesian_to_equinoctial(out.X(1:3,n),out.X(4:6,n));
        lMmean = mod(lMmean,twopi);    
        out.E(:,n) = [nmean,afmean,agmean,chimean,psimean,lMmean]';

        % Calculate Jacobian matrix and inverse
        out.J(:,:,n) = ...
            jacobian_equinoctial_to_cartesian(out.E(:,n),out.X(:,n));
        out.K(:,:,n) = out.J(:,:,n)\I6x6;

        % Calculate equinoctial covariance, enforcing symmetry to correct
        % round-off error effects
        out.Q(:,:,n) = cov_make_symmetric( ...
            out.K(:,:,n) * out.P(:,:,n) * out.K(:,:,n)' );

    end
    
end

return
end

% =========================================================================

function [T,Nexcl] = combine_eph_times(T1,T2,StepFracMin,verbose)

% Combine two sets of ephemeris table times. Specifically, combine 
% the set in array T2 with the set in array T1, excluding those within a
% small fraction of an time step of the original set T1

% Number of time points in second set
N2 = numel(T2); % NminTimes = numel(T1);

excl = false(size(T2)); Nexcl = 0;

idx = (T2(1) <= T1) & (T1 <= T2(N2));
idx = find(idx); Nidx = numel(idx);

for ii = 1:Nidx
    
    i = idx(ii);
    
    [c,b,coincidence] = find_bracketing_indices(T1(i),T2);

    if coincidence
        % Minimum perfectly coincides with an ephemeris point
        excl(c) = true; Nexcl = Nexcl+1;
    else
        StepFrac = (T1(i)-T2(c))/(T2(b)-T2(c));
        if StepFrac < 0 || StepFrac > 1
            error('Invalid time-step fraction calculated');
        end
        if StepFrac < StepFracMin
            excl(c) = true;  Nexcl = Nexcl+1;
        end
    end
        
end

if verbose
    disp([' Number of eph2 times excluded as too close to eph1 times: ' ...
        num2str(Nexcl)]);
end

% Combine the ephemeris tables, excluding eph points as required
ndx = ~excl;    
T  = cat( 2 , T1 , T2(ndx) );

% Sort the combined ephemeris
T = sort(T);

return
end

% =========================================================================

function integrand = Ncdot_integrand(rht,mur,muv,Ainv,R,Q1,Q2,logZ,slow_method)

% Calculate the integrand for the Ncdot unit-sphere integral.

% Initialize output

sz = size(rht);
Nsz = numel(sz);

if (Nsz == 2)
    sz(1) = 1;
else
    sz = sz(2:end);
end

integrand = NaN(sz);

% Slow vs fast method

if slow_method

    % Number of elements in (ph,u) arrays

    N = prod(sz);

    % Check if velocity dispersion will be zero. Here 
    %   Q2 = C - B * Ainv * B'.
    
    zero_sig2 = (max(abs(Q2(:))) == 0);

    % Perform calculate for each of the elements
    
    sqrt2 = sqrt(2);
    sqrtpi = sqrt(pi);
    sqrt2pi = sqrt2*sqrtpi;
    
    num_nonpos_sig2 = 0;

    for n=1:N

        % Calculate rhat

        rhat = rht(:,n);
        
        % Calculate difference from center
        
        dr = R*rhat-mur;
        
        % Calculate the nu(rhat,t) factor.
        
        if zero_sig2
            
            % Limiting nu factor  as sigma -> 0+ (C12a eq 39)
            
            nu = max(0, -muv' * rhat);
            
        else

            % Positive sigma C12a eqs 31, 36 & 37. Here
            %   Q1 = B * Ainv;
            %   Q2 = C - B * Ainv * B'
        
            sig2 = rhat' * Q2 * rhat;

            % Handle nonpositive and positive sig2 values
            
            if (sig2 <= 0)
                % Limiting nu factor  as sigma -> 0+ (C12a eq 39)
                num_nonpos_sig2 = num_nonpos_sig2+1;
                nu = max(0, -muv' * rhat);
            else
                sig = sqrt(sig2);
                nu0 = rhat' * (muv + Q1 * dr);
                nus = nu0 / sig / sqrt2;
                H = exp(-nus^2)-sqrtpi*nus*erfc(nus);
                nu = sig * H / sqrt2pi;
            end

        end
        
        % Calculate the Mahalanobis distance squared

        MD2 = dr' * Ainv * dr;

        % Calculate the integrand, which is the MVN N3 function 

        neglogN3 = logZ + 0.5*MD2;
        integrand(n) = nu * exp(-neglogN3);

        % disp(['n, ph, u, i = ' num2str(n) ' ' num2str(ph(n)) ...
        %       ' ' num2str(u(n)) ' ' num2str(integrand(n))]);
          
        if isnan(integrand(n)) || isinf(integrand(n)) || (integrand(n) < 0)
            error('Bad Ncdot integrand');
        end
        
    end
    
    if (num_nonpos_sig2 > 0)
        warning([num2str(num_nonpos_sig2) ' of ' ...
            num2str(N) ...
            ' velocity sigma^2 factors are not positive.']);
    end
    
else
    
    % Allocate space ultimately for difference-from-center vectors,
    % i.e. dr = r-mur.  NOTE: initially this will be used to hold rhat
    % vectors, and later changed to dr vectors.
    
    sznew = cat(2,3,cat(2,1,sz));
    rhat = zeros(sznew);
    
    szmat = cat(2,1,cat(2,1,sz));    
    
    % Expand the rhat vector array
    
    rhat(1,1,:,:) = rht(1,:,:);
    rhat(2,1,:,:) = rht(2,:,:);
    rhat(3,1,:,:) = rht(3,:,:);
    
    % Calculate the dr = r-mur vectors
    
    dr = R*rhat - repmat(mur,szmat);
    
    % Check if velocity dispersion will be zero. Here 
    %   Q2 = C - B * Ainv * B'.
    
    zero_sig2 = (max(abs(Q2(:))) == 0);
    
    % Calculate the nu(rhat,t) factor.

    if zero_sig2

        % Limiting nu factor  as sigma -> 0+ (C12a eq 39)

        nu = -squeeze(multiprod(multitransp(repmat(muv,szmat)),rhat));
        
        % Ensure nu is not negative

        nu(nu < 0) = 0;

    else
        
        % Use the rhat vectors to calculate the nu(rhat,t) factors, 
        % using C12a eqs 31, 36 and 37.  Here
        %
        %   Q1 = B * Ainv;
        %   Q2 = C - B * Ainv * B'.

        sig2 = squeeze(multiprod( ...
            multitransp(rhat),multiprod(repmat(Q2,szmat),rhat)));
        
        % Calculate nu values.  First find any nonpositive sigma^2 values,
        % and mark the associated sigma values with NaNs.

        nonpos_sig2 = (sig2(:) <= 0);        
        num_nonpos_sig2 = sum(nonpos_sig2);
        any_nonpos_sig2 = (num_nonpos_sig2 > 0);
        
        if any_nonpos_sig2
            sig = NaN(size(sig2));
            pos_sig2 = ~nonpos_sig2;
            sig(pos_sig2) = sqrt(sig2(pos_sig2));
            % warning([num2str(num_nonpos_sig2) ' of ' ...
            %     num2str(numel(sig2(:))) ...
            %     ' velocity sigma^2 factors are not positive.']);
        else
            sig = sqrt(sig2);
        end
        
        nu0 = squeeze(multiprod(multitransp(rhat),repmat(muv,szmat) ...
            + multiprod(repmat(Q1,szmat),dr)));
        
        % Use nu0 to calculate nu
        
        sqrt2 = sqrt(2);
        sqrtpi = sqrt(pi);
        sqrt2pi = sqrt2*sqrtpi;
        
        nus = nu0 ./ (sig * sqrt2);
        H = exp(-nus.^2) - sqrtpi * (nus .* erfc(nus) );
        nu = sig .* H / sqrt2pi;

        % Account for the nonpositive sigma^2 values
        
        if any_nonpos_sig2
            
            if (Nsz == 2)
                
                nu_nonpos_sig2  = -muv' * rht(:,nonpos_sig2);
                nu_nonpos_sig2(nu_nonpos_sig2 < 0) = 0;
                nu(nonpos_sig2) = nu_nonpos_sig2;
                
            else
                
                % Size of nu array

                sznu = size(nu);

                % Find the single-index indices for the nonpositive sigma^2
                % values

                ndx_nonpos_sig2 = find(nonpos_sig2);

                % Loop over indices and fix the bad nu values resulting from
                % nonpositive sigma^2 values

                for nn=1:num_nonpos_sig2

                    % Get subscripts mapping into nu array

                    kk = ndx_nonpos_sig2(nn);
                    [ii,jj] = ind2sub(sznu,kk);

                    % Limiting nu factor  as sigma -> 0+ (C12a eq 39)                

                    rh = rht(:,ii,jj);

                    nu(ii,jj) = max(0, -muv' * rh);

                end
                
            end

        end
        
    end
    
    % Calculate the Mahalanobis distance squared
    
    MD2 = squeeze(multiprod(multitransp(dr),multiprod(repmat(Ainv,szmat),dr)));
    
    % Calculate the integrand, which is the product of the MVN N3 function
    % and the nur(r,t) factor (C12a eq 36a).
    
    neglogN3 = logZ + 0.5*MD2;
    integrand = nu .* exp(-neglogN3);
    
end

return
end

% =========================================================================

function integrand = Ncdot_quad2d_integrand(ph,u,mur,muv,Ainv,R,Q1,Q2,logZ,slow_method)

% Calculate the integrand for the Ncdot unit-sphere integral, assuming that 
% (ph,u) are 2D matrices of equal dimension, as would be required for use
% with Matlab's quad2d function.

% Get the size of the (ph,u) arrays

sz = size(ph);

% Allocate space for rhat vectors

sznew = cat(2,3,sz);
rhat = zeros(sznew);

% Complement of u

up = real(sqrt(1-u.^2));

% Calculte r0-hat vectors

rhat(1,:,:) = cos(ph) .* up;
rhat(2,:,:) = sin(ph) .* up;
rhat(3,:,:) = u;

% Calculate integrand

integrand = Ncdot_integrand(rhat,mur,muv,Ainv,R,Q1,Q2,logZ,slow_method);

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
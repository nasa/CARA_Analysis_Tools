function [PcOne,Diluted,PcMax,SfMax,Pc,Sf,conv,iter] = PcDilution(r1,v1,cov1,r2,v2,cov2,HBR,params)
% PcDilution - Evaluates collision probability (Pc) dilution and calculates
%              the associated maximum Pc value, resulting from scaling a
%              conjunction's secondary (default) or primary 3x3 position
%              uncertainty covariance matrix by applying a "sigma scaling
%              factor" (Sf). 
%
% Syntax: [Diluted,PcMax,SfMax,PcOne,Pc,Sf,conv,iter] = PcDilution(r1,v1,cov1,r2,v2,cov2,HBR,params)
%
% Inputs:
%    r1      - Primary object's position vector in ECI coordinates
%              [1x3 or 3x1] (meters)
%    v1      - Primary object's velocity vector in ECI coordinates
%              [1x3 or 3x1] (meters/second)
%    cov1    - Primary object's covariance matrix in ECI coordinate frame
%              (3x3 or 6x6) (meters and meters/s)
%    r2      - Secondary object's position vector in ECI coordinates
%              [1x3 or 3x1] (meters)
%    v2      - Secondary object's velocity vector in ECI coordinates
%              [1x3 or 3x1] (meters/second)
%    cov2    - Secondary object's covariance matrix in ECI coordinate frame
%              (3x3 or 6x6) (meters and meters/s)
%    HBR     - Hard body region (meters)
%    params  - Run parameters [optional]
%            - params.RefineCA = Flag to refine CA point before analysis
%               (default = false; assumes true CA quantities are input)
%            - params.PriSecScaling = Scaled scovariance, can be 'primary'
%               'secondary' or 'both' (default = 'secondary')
%            - params.RedFact = Pc reduction factor required for
%               convergence. Recommendations:
%               1) Use RedFact = 10 to 100 when plotting the (Sf,Pc)
%                  output curve.
%               2) Use RedFact = 1.1 (the default) for automated processing
%                  with no post facto plotting.
%            - params.ConstTol - Tolerance to evaluate the occurence of
%                constant Pc values as a function of the scale factor Sf
%                (default = 1e-6).
%            - params.DilutTol - Tolerance to evaluate the occurence of
%                a diluted PcMax value (default = 1e-3).
%            - params.itermax - Maximum number of iterations allowed to
%               converge upon PcMax (default = 100)
%            - params.SfInit - Initial range of sigma-value scale factors
%               (default = [0.5 2.0]).
%            - params.dLfLevel - Log10 sigma scale factor spacing, initial
%               and refinement levels (default = 2*[1e-2 1e-3 1e-4]).
%            - params.verbose - Flag for verbose operation, used for
%               development and debugging (default = false).
%            - params.plotting - Plotting monitorint level, used for
%               development and debugging:
%               0 = No plotting (default)
%               1 = Plot final results only
%               2 = Plot iteration-by-iteration results
%
% Outputs:
%    PcOne   - Pc calculated for unscaled covariances (i.e., for a sigma
%              scale factor of one, or Sf = 1), commonly referred to
%              as "the conjunction's nominal 2D-Pc value."
%    Diluted - Flag to indicate probability is in dilution region.
%    PcMax   - Maximum Pc value found for any value of applied sigma scale
%              factor (Sf).
%    SfMax   - Sigma scale factor for which maximum Pc occurs.
%    Pc      - The Pc buffer created in the iterative search for PcMax.
%    Sf      - Sigma scale factor buffer created in the search for PcMax.
%              (Plotting y=Pc vs x=Sf creates the Pc dilution curve.)
%    conv    - Flag indicating satisfactory convergence to find PcMax.
%    iter    - Number of iterations performed in finding PcMax.
%
% References:  M.Hejduk (2019) "Satellite Conjunction Assessment Risk 
%              Analysis for 'Dilution Region' Events: Issues and
%              Operational Approaches" Space Trafic Managment Conference,
%              28, https://commons.erau.edu/stm/2019/presentations/28
%
%              S.Alfano (2005) "Relating Position Uncertainty to Maximum
%              Conjunction Probability" The Journal of the Astronautical
%              Sciences, Vol.53, No.2, pp.193-205.
%
% Example/Validation Cases:
%
% Other m-files required:
%   FindNearbyCA.m
%   isconstant.m
%   PcElrod.m
% Subfunctions: None
% MAT-files required: None
% See also: None
%
% November 2019; Last revision: 2019-Nov-12
%
% ----------------- BEGIN CODE -----------------

    % Set up default parameters
    
    Nargin = nargin;

    if Nargin < 8
        params = [];
    end
    
    % Refine CA point before analysis
    
    if ~isfield(params,'RefineCA') || isempty(params.RefineCA)
        params.RefineCA = false;
    end
    
    % Scaling secondary covariance (default), primary or both
    
    if ~isfield(params,'PriSecScaling') || isempty(params.PriSecScaling)
        params.PriSecScaling = 'Secondary';
    end
    
    % Pc reduction factor required for convergence
    
    if ~isfield(params,'RedFact') || isempty(params.RedFact)
        params.RedFact = 1.1;
    end

    % PcMax constance tolerance
    
    if ~isfield(params,'ConstTol') || isempty(params.ConstTol)
        params.ConstTol = 1e-6;
    end
    
    % PcMax dilution tolerance
    
    if ~isfield(params,'DilutTol') || isempty(params.DilutTol)
        params.DilutTol = 1e-3;
    end
    
    % Maximum iteration numbers require for convergence
    
    if ~isfield(params,'itermax') || isempty(params.itermax)
        params.itermax = 100;
    end

    % Initial sigma scale factors and log10 values
    
    if ~isfield(params,'SfInit') || isempty(params.SfInit)
        params.SfInit = [0.5 2.0];
    end
    
    % Log10 sigma scale factor spacing, initial and refinement levels
    
    if ~isfield(params,'dLf') || isempty(params.dLf)
        params.dLfLevel = 2*[1e-2 1e-3 1e-4];
    end
    
    % Printing and plotting (mostly for debugging)
    
    if ~isfield(params,'verbose') || isempty(params.verbose)
        params.verbose = false;
    end
    
    if ~isfield(params,'plotting') || isempty(params.plotting)
        params.plotting = 0;
    end
    
    % Check for valid parameters
    
    SfInit = unique(params.SfInit); % Sort into increasing order
    if (numel(SfInit) <= 1)
        error('Invalid SfInit parameter; [0.5 2.0] recommended');
    end
    
    RedFact = params.RedFact;
    if (RedFact <= 1)
        error('Invalid RedFact parameter; must be > 1');
    end
    
    PriSecScaling = params.PriSecScaling;
    if ischar(PriSecScaling)
        PriSecScaling = lower(PriSecScaling);
        switch PriSecScaling
            case {'primary','pri'}
                PriSecScaling = [true false];
            case {'secondary','sec'}
                PriSecScaling = [false true];
            case {'both'}
                PriSecScaling = [true true];
            otherwise
                error('Invalid PriSecScaling string parameter');
        end
    else
        if ~islogical(PriSecScaling)   || ...
           (numel(PriSecScaling) ~= 2) || ...
           ~any(PriSecScaling)
            error('Invalid PriSecScaling logical parameter');
        end
    end
    
    % Reshape the input vectors to be 1x3

    r1 = reshape(r1,1,3); v1 = reshape(v1,1,3);
    r2 = reshape(r2,1,3); v2 = reshape(v2,1,3);
    
    % Refine the CA point, if required
    
    if params.RefineCA
        [~,X1,X2] = FindNearbyCA([r1 v1]',[r2 v2]');
        r1 = X1(1:3)'; v1 = X1(4:6)';
        r2 = X2(1:3)'; v2 = X2(4:6)';
    end
    
    % Ensure the covariance matrices are 3x3 (i.e., trim any 6x6 matrices)
    
    cov1 = cov1(1:3,1:3); cov2 = cov2(1:3,1:3);

    % Initialize output buffers to null, just in case no iterations are
    % required to converge upon dilution solution
    
    Pc = [];  Sf = []; conv = true; iter = 0;
    
    % If both covariances are zero then there can be no Pc dilution, and
    % the Pc estimate reduces to that of a hard-body collision
    
    iszerocov1 = all(cov1(:) == 0); iszerocov2 = all(cov2(:) == 0);
    
    if iszerocov1 && iszerocov2
        if params.verbose
            disp('Both covariances zero => no dilution possible');
        end
        Diluted = false;
        DCA = norm(r1-r2); % Distance of closest approach
        if (DCA < HBR)
            PcOne = 1; % Hard-body collision
        else
            PcOne = 0; % No hard-body collision
        end
        SfMax = 1; PcMax = PcOne;
        return;
    end
    
    % Find the Pc value for no sigma scaling (i.e., for Sf = 1)
    
    [PcOne,CpOne] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
    
    % Return if 2D-Pc is NaN, but indicate no convergence
    
    if isnan(PcOne)
        conv = false; Diluted = NaN; SfMax = NaN; PcMax = NaN;
        return;
    end
    
    % Compare scaling factor to analytical approximations
    if (params.plotting > 0) && all(PriSecScaling)
        [Evec,Sig2] = eig(CpOne);
        Sig2 = diag(Sig2);
        [MaxSig2,nmx] = max(Sig2);
        [MinSig2,~] = min(Sig2);
        theta = atan2(Evec(2,nmx),Evec(1,nmx));
        cth = cos(theta); sth = sin(theta);
        dist = norm(r1-r2);
        sigx = sqrt(MinSig2);
        sigy = sqrt(MaxSig2);
        smallHBR = HBR < 0.2*min([dist,sigx,sigy]);
        if smallHBR
            SfMaxApp = (dist/(sqrt(2)*sigx*sigy))*sqrt((sigx*cth)^2+(sigy*sth)^2);
            PcOneApp = (HBR^2/2/sigx/sigy)*exp(-0.5*dist^2*((sth/sigx)^2+(cth/sigy)^2));
            PcMaxApp = (HBR^2/2/sigx/sigy/SfMaxApp^2)*exp(-0.5*dist^2*((sth/SfMaxApp/sigx)^2+(cth/SfMaxApp/sigy)^2));
            if params.verbose
                disp(['Approximation for HBR << min(MissDist,Sigx,Sigy) or ' ...
                    num2str(HBR) ' << ' num2str(min([dist,sigx,sigy]))]);
                disp(['SfMaxApp = ' num2str(SfMaxApp) ...
                     ' PcMaxApp = ' num2str(PcMaxApp) ...
                     ' Pc2DApp = ' num2str(PcOneApp)]);
            end
        end
    end
    
    % Check if the only scaled covariance is zero
    
    if iszerocov1 && PriSecScaling(1) && ~PriSecScaling(2)
        % No dilution because zero primary cov. can't be scaled down
        if params.verbose
            disp('Zero primary covariance => no primary-scaling dilution possible');
        end
        Diluted = false; SfMax = 1; PcMax = PcOne;
        return;
    elseif iszerocov2 && PriSecScaling(2) && ~PriSecScaling(1)
        % No dilution because zero secondary cov. can't be scaled down
        if params.verbose
            disp('Zero secondary covariance => no secondary-scaling dilution possible');
        end
        Diluted = false; SfMax = 1; PcMax = PcOne;
        return;
    end
    
    % Find the Pc value for complete sigma down-scaling (i.e., Sf = 0)

    if all(PriSecScaling)
        % Scaling both pri & sec covariance to zero again reduces to a
        % hard-body collision probability in the small scale factor limit,
        % because both covariances approach zero
        DCA = norm(r1-r2); % Distance of closest approach
        if (DCA < HBR)
            PcZero = 1; % Hard-body collision
        else
            PcZero = 0; % No hard-body collision
        end
    else
        % Scale pri or sec covariance using Sf = 0
        if PriSecScaling(1); f1 = 0; else; f1 = 1; end
        if PriSecScaling(2); f2 = 0; else; f2 = 1; end
        % Pc value for Sf = 0
        PcZero = PcElrod(r1,v1,f1*cov1,r2,v2,f2*cov2,HBR);
    end
    
    % Initial sigma scale factors and log10 values
    
    Sf1 = params.SfInit(1);
    Sf2 = params.SfInit(2);
    
    Lf1 = log10(Sf1);
    Lf2 = log10(Sf2);
    
    log10_two = log10(2);
    
    % Log10 sigma scale factor spacing

    NLfLevel = numel(params.dLfLevel);
    nLfLevel = 1;
    dLf = min(params.dLfLevel(nLfLevel),(Lf2-Lf1)/5);
    
    % Initialize buffers
    
    % Initialize plotting
    
    if params.plotting > 0
        figure;
    end
    
    % Iterate until converged or iteration limit exceeded
    
    iterating = true; conv = false; Pc0 = PcOne;

    while iterating
        
        % Log10 sigma scale factor spacing
    
        Lf = (Lf1:dLf:Lf2)';
        N = numel(Lf);

        % Scaled covariances

        Sfb = 10.^linspace(Lf1,Lf2,N)'; f = Sfb.^2;
        if PriSecScaling(1); f1 = f; else; f1 = ones(size(f)); end
        if PriSecScaling(2); f2 = f; else; f2 = ones(size(f)); end

        c1 = NaN(3,3,N); c2 = c1;
        
        for n=1:N
            c1(:,:,n) = cov1*f1(n);
            c2(:,:,n) = cov2*f2(n);
        end

        % Pc for scaled covariances

        Nx1 = [N 1];
        Pcb = PcElrod( ...
            repmat(r1,Nx1),repmat(v1,Nx1),c1, ...
            repmat(r2,Nx1),repmat(v2,Nx1),c2, ...
            repmat(HBR,Nx1));
        
        % Add to buffers
        
        Sf = cat(1,Sf,Sfb); Pc = cat(1,Pc,Pcb); Nf = numel(Sf);
        
        % Sort in increasing scale factor
        
        if iter > 0
            [Sf,srt] = sort(Sf); Pc = Pc(srt);
        end
        
        % Check if scale factor extremes are tolerable
        
        Pc0 = max(Pc);
        PcRed = Pc0/RedFact;
        
        if (Pc0 == 0)
            % Expand both bounds because all scale factors explored so far
            % have resulted in Pc = 0
            if params.verbose
                disp('Expanding both scale factor bounds');
            end
            Lf1 = Lf1-log10_two;
            Lf2 = Lf2+log10_two;
        elseif (Pc(1) > PcRed) && ~isconstant([Pc(1) PcZero],params.ConstTol)
            % Decrease lower bound
            if params.verbose
                disp('Decreasing lower scale factor');
            end
            Lf2 = log10(Sf(1))-dLf;
            Lf1 = Lf2-log10_two;
        elseif (Pc(end) > PcRed) && ~isconstant(Pc(Nf-N+1:Nf),params.ConstTol)
            % Increase upper bound
            if params.verbose
                disp('Increasing upper scale factor');
            end
            Lf1 = log10(Sf(end))+dLf;
            Lf2 = Lf1+log10_two;
        elseif nLfLevel < NLfLevel
            % Fine tuning maximum Pc
            nLfLevel = nLfLevel+1;
            if params.verbose
                disp(['Fine tuning level ' num2str(nLfLevel)]);
            end
            [~,imax] = max(Pc);
            i1 = max(imax-1,1);
            i2 = min(imax+1,numel(Pc));
            Lf1 = log10(Sf(i1));
            Lf2 = log10(Sf(i2));
            dLf = min(params.dLfLevel(nLfLevel),(Lf2-Lf1)/5);
        else
            % Converged
            if params.verbose
                disp('Converged');
            end
            iterating = false; conv = true;
        end
        
        % Increment iteration counter and terminate if too large
        
        iter = iter+1;
        
        if (iter > params.itermax)
            iterating = false;
        end
        
        if params.plotting > 1
            clf;
            loglog(Sf,Pc,'+--k');
            hold on;
            plot(1,PcOne,'or');
            hold off;
            drawnow;
            keyboard;
        end
        
    end
    
    % Add unscaled Pc to output buffers, and make unique/sorted 
    
    Sf = cat(1,1,Sf);
    Pc = cat(1,PcOne,Pc);
    [Sf,srt] = unique(Sf); Pc = Pc(srt);

    % Define best-estimate Pc maximum

    [PcMax,imax] = max(Pc);
    if isnan(PcMax)
        SfMax = NaN;
        Diluted = NaN;
    else
        SfMax = Sf(imax);
        Diluted = SfMax < 1; % Dilution occurs only for sigma scaling < 1
        if Diluted
            % Dilution flagged only when PcMax/Pc-2D exceeds the 
            % specified dilution tolerance
            Diluted = ~isconstant([PcOne PcMax],params.DilutTol);
        end
        if params.verbose
            disp(['Primary scaling=' num2str(PriSecScaling(1)) ...
               ' Secondary scaling=' num2str(PriSecScaling(2)) ...
               ' Diluted=' num2str(Diluted) ...
               ' SfMax=' num2str(SfMax) ...
               ' PcMax=' num2str(PcMax) ...
               ' Pc2D=' num2str(PcOne)]);
        end
    end
    
    % Issue failure to convergence warning
    
    if params.verbose && ~conv
        warning(['PcDilution failed to converge: iter = ' num2str(iter)]);
    end
    
    if params.plotting > 0
        if all(PriSecScaling) && smallHBR
            smallHBRapp = true;
            lsty = '--';
        else
            smallHBRapp = false;
            lsty = '-';
        end
        clf;
        loglog(Sf,Pc,[lsty 'k']);
        hold on;
        plot(1,PcOne,'or');
        plot(SfMax,PcMax,'sm');
        hold off;
        if smallHBRapp
            ssx=Sf*sigx; ssy=Sf*sigy;
            PcApp = ((HBR^2)./2./ssx./ssy).*exp(-0.5*dist^2*((sth./ssx).^2+(cth./ssy).^2));
            hold on;
            plot(Sf,PcApp,':b');
            hold off;
            yrng = plot_range([min(Pc(Pc > 0)) max(Pc) PcMax],0);
            ylim(yrng);
        end
        xlabel('Sigma Scale Factor');
        ylabel('Pc');
        title(['Primary scaling = ' num2str(PriSecScaling(1)) ...
            ' Secondary scaling = ' num2str(PriSecScaling(2))]);
        keyboard;
    end
    
    return;
    
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D.Hall         | 2019-SEP-04 | Initial Development
% D.Hall         | 2019-SEP-09 | Added DilutTol parameter
% D.Hall         | 2019-NOV-12 | Added code to handle cases where all Pc
%                                values are zero within the initial range
%                                of sigma-value scale factors.
%

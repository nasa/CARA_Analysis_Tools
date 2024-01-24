function [PcXC,CovXC,DCPvalues] = Pc_CovXCorr(cdmhead,cdmobj,HBR,DCPoption,verbose)
% Pc_CovXCorr - Calculate a 2D-Pc estimate with a joint relative position
%               covariance matrix corrected for cross correlation effects,
%               using the algorithm of Casali et al (2018) AAS 18-272.
%
% Syntax: [PcXC,CovXC,DCPvalues] = Pc_CovXCorr(cdmhead,cdmobj,HBR,DCPoption,verbose)
%
% Inputs:
%
%    cdmhead    - CDM header structure as returned by function "read_cdm.m"
%                 [1x1 struct]
%    cdmobj     - CDM object structures as returned by "read_cdm.m"
%                 [2X1 struct]
%    HBR        - Hard body radius for 2D-Pc calculation (m)
%                 [1x1 double]
%    DCPoption  - Integer specifying how to get required DCP values when
%                 applying the Casali et al (2018) algorithm (default = 1):
%                   1 - Use DCP Density Forecast Uncertainty and 
%                       DCP Sensitivity Vector values, but only if they are
%                       specified in the CDM explicitly. (These DCP values
%                       are written into CDMs created by ASW version 19.2
%                       and after.)
%                   2 - Use DCP Density Forecast Uncertainty and 
%                       DCP Sensitivity Vector values if they are
%                       specified in the CDM, but otherwise use the
%                       EDR approximation for these DCP values.
%                   3 - Use DCP Density Forecast Uncertainty and 
%                       DCP Sensitivity Vector value EDR approximations,
%                       whether the actual DCP are specified in the CDM
%                       or not.
%                 WARNING: USE OF DCPoption = 1 (THE DEFAULT) IS THE ONLY
%                          MODE RECOMMENDED FOR RISK ANALYSIS OF ACTUAL
%                          CONJUNCTIONS. USING DCPoption = 2 OR 3 CAN 
%                          POTENTIALLY EMPLOY DCP VALUES ESTIMATED USING
%                          A VERY ROUGH EDR APPROXIMATION, SOMETIMES 
%                          RESULTING IN NPD CovXC MATRICES AND INACCURATE
%                          PcXC RESULTS. USING DCPoption = 2 OR 3 IS ONLY
%                          RECOMMENDED FOR DEVELOPMENT AND DEBUGGING.
%   verbose     - Verbose operation (default = false)
%
% Outputs:
%
%   PcXC        - The 2D-Pc estimate obtained by applying the
%                 Casali et al (2018) algorithm to correct for covariance
%                 cross correlation effects. Note, if no DCP values were
%                 available during the processing (either as specified in
%                 CDM or approximated if permitted), then PcXC is returned
%                 as an empty set [].
%   CovXC       - The 3x3 relative position covariance obtained using the
%                 Casali et al (2018) algorithm to correct for covariance
%                 cross correlation effects. Note, if no DCP values were
%                 available during the processing (either as specified in
%                 CDM or approximated if permitted), then CovXC is returned
%                 as an empty set [].
%   DCPvalues   - A structure holding the DCPvalues used in the calculation
%                 for both the primary and the secondary objects.
%
% References:
%
%    S.J.Casali et al (2018) "Effect of Cross-Correlation of Orbital Error
%    on Probability of Collision Determination" AAS 18-272.
%    (Hereafter refered to as C18.)
%
% Examples/Validation Cases: TBD
%
% Other m-files required:
%   cdmobj_J2Kstatecov.m
%   cdmobj_DCPvalues.m
%   RelAtmoDensUnc.m
%   convert_cartesian_to_equinoctial.m
%
% Subfunctions: None
%
% MAT-files required: None
%
% See also: None
%
% September 2019; Last revision: 2019-Sep-30
%
% ----------------- BEGIN CODE -----------------

    % Initializations and defaults

    % Initialize input
    Nargin = nargin;
    if (Nargin < 5) || isempty(verbose)
        verbose = false;
    end
    if (Nargin < 4) || isempty(DCPoption)
        DCPoption = 1;
    end
    
    % Perform covariance cross correlation processing

    % Get the J2K states and covariances from the CDM data
    [Xp,Pp] = cdmobj_J2Kstatecov(cdmhead,cdmobj(1));
    [Xs,Ps] = cdmobj_J2Kstatecov(cdmhead,cdmobj(2));

    % Extract the DCP values for the primary and secondary from the CDM
    % data, if required.  These are the sigma values and G vectors from
    % C18 eq (11).
    if DCPoption ~= 3
        % Extract the DCP values from the CDM object data structures
        [sigp,Gvecp] = cdmobj_DCPvalues(Xp,cdmobj(1));
        [sigs,Gvecs] = cdmobj_DCPvalues(Xs,cdmobj(2));
        % If the primary and secondary have an incomplete set of DCP
        % values, then issue a warning and ignore both sets
        isemptyDCPp = isempty(sigp) | isempty(Gvecp);
        isemptyDCPs = isempty(sigs) | isempty(Gvecs);
        if ~isequal(isemptyDCPp,isemptyDCPs)
            if verbose
                disp(['Incomplete set(s) of primary or secondary ' ...
                      'DCP values found in CDM object data; ' ...
                      'neglecting both sets.']);
            end
            isemptyDCPp = true;
            isemptyDCPs = true;
        end
        % % Temporary code to handle wrong units for G vectors
        % if ~isemptyDCPp
        %     Gvecp = 1e3*Gvecp;
        %     warning('Multiplying Gvecp by 1000 to convert from km to m units');
        % end
        % if ~isemptyDCPs
        %     Gvecs = 1e3*Gvecs;
        %     warning('Multiplying Gvecs by 1000 to convert from km to m units');
        % end
    else
        % DCPoption = 3 specifies not to use DCP values whether they exist
        % in the CDM data or not, so keep them empty at this point
        isemptyDCPp = true;
        isemptyDCPs = true;
    end
   
    % Initialize the approximation flag
    DCPapprox = false;

    % If the primary and secondary DCP values are nonexistent (i.e.,
    % empty), then estimate them using the EDR-based approximation given
    % by C18 in eqs (17)-(20), but only if permitted by the DCPoption
    if isemptyDCPp && isemptyDCPs && (DCPoption ~= 1)
        
        % Ensure that the TIME_LASTOB_END values exist
        if isfield(cdmobj(1),'TIME_LASTOB_END') && ...
          ~isempty(cdmobj(1).TIME_LASTOB_END)   && ...
           isfield(cdmobj(2),'TIME_LASTOB_END') && ...
          ~isempty(cdmobj(2).TIME_LASTOB_END)
        
            % Set the approximation flag
            DCPapprox = true;

            % Earth gravitational constant (EGM-96) [m^3/s^2]
            mu  = 3.986004418e14;

            % Earth equatorial radius (m)
            Re = 6378.137e3;        

            % Extract TCA from the CDM header
            TCA = datenum(cdmhead.TCA, 'yyyy-mm-ddTHH:MM:SS.FFF');

            % Extract CDM last observation epochs
            TEp = datenum(cdmobj(1).TIME_LASTOB_END, 'yyyy-mm-ddTHH:MM:SS.FFF');
            TEs = datenum(cdmobj(2).TIME_LASTOB_END, 'yyyy-mm-ddTHH:MM:SS.FFF');

            % Approximate the DCP values sig and Gvec for primary
            t = (TCA-TEp)*86400;
            EDR = cdmobj(1).SEDR;
            % EDR = max(0,EDR);
            rvec = Xp(1:3)';
            vvec = Xp(4:6)';
            alt_m = norm(rvec)-Re;
            sigp = RelAtmoDensUnc(alt_m/1e3);
            [a,n] = convert_cartesian_to_equinoctial(rvec,vvec,1,mu);
            DM = 1.5*EDR*(t^2)/sqrt(mu*a);
            Gvecp = NaN(1,6);
            Gvecp(1:3) = (DM/n)*vvec';
            isemptyDCPp = false;

            % Approximate the DCP values sig and Gvec for secondary
            t = (TCA-TEs)*86400;
            EDR = cdmobj(2).SEDR;
            % EDR = max(0,EDR);
            rvec = Xs(1:3)';
            vvec = Xs(4:6)';
            alt_m = norm(rvec)-Re;
            sigs = RelAtmoDensUnc(alt_m/1e3);
            [a,n] = convert_cartesian_to_equinoctial(rvec,vvec,1,mu);
            DM = 1.5*EDR*(t^2)/sqrt(mu*a);
            Gvecs = NaN(1,6);
            Gvecs(1:3) = (DM/n)*vvec';
            isemptyDCPs = false;
            
        end
     
    end
    
    % If the DCPvalues remain empty at this point, then covariance cross 
    % correlation effects cannot be calculated, so return empty-set output
    if isemptyDCPp && isemptyDCPs
        if verbose
            disp(['No DCP values available, ' ...
                  'so no cov. cross correlation correction performed'])
        end
        PcXC      = [];
        CovXC     = [];
        DCPvalues = [];
        return;
    end
    
    % Define the output DCP values
    DCPvalues.sigp = sigp; DCPvalues.Gvecp = Gvecp;
    DCPvalues.sigs = sigs; DCPvalues.Gvecs = Gvecs;
    DCPvalues.approximated = DCPapprox;
    
    % Warn user about EDR approximation
    if DCPvalues.approximated && DCPoption ~= 3
        warning(['DCP values estimated using the EDR approximation, ' ...
                 'potentially yielding very inaccurate results']);
    end

    % Extract the position part of the sensitivity vectors,
    % and make column vectors
    Gp = Gvecp(1:3)'; 
    Gs = Gvecs(1:3)';

    % Calculate the relative position covariance corrected for covariance
    % cross correlation, CovXC (denoted as Pm in C18 eq. 11)
    CovXC = Ps(1:3,1:3) + Pp(1:3,1:3) - (sigs*sigp) * (Gs*Gp'+Gp*Gs');

    % Calculate the 2D-Pc value using the corrected covariance
    PcXC = PcElrod(Xp(1:3),Xp(4:6),CovXC,      ...
                   Xs(1:3),Xs(4:6),zeros(3,3), ...
                   HBR);

    % Report results if required
    if verbose
        % Calculate the 2D-Pc neglecting any covariance correlation
        Pc = PcElrod(Xp(1:3),Xp(4:6),Pp(1:3,1:3), ...
                     Xs(1:3),Xs(4:6),Ps(1:3,1:3), ...
                     HBR);
        % Report Pc results
        PcFmt = '%0.7e';
        disp(['Pc neglecting cov. cross corr. = ' num2str(Pc  ,PcFmt) ]);
        disp(['Pc correcting cov. cross corr. = ' num2str(PcXC,PcFmt) ]);
        if Pc == PcXC
            rat = 1;
        else
            rat = PcXC/Pc;
        end
        disp(['   corrected/uncorrected ratio = ' num2str(rat,PcFmt)]);
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
% D.Hall         | 2019-SEP-30 | Initial Development
%

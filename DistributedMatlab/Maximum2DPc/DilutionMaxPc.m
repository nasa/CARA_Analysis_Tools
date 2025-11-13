function [PcMax,Diluted,Pri,Sec] = DilutionMaxPc(r1,v1,cov1,r2,v2,cov2,HBR,params)
% DilutionMaxPc - Evaluates collision probability (Pc) dilution and
%                 the associated maximum Pc value resulting from scaling a
%                 conjunction's 3x3 secondary or 3x3 primary covariance, or
%                 combined covariance (default).
%
% Syntax: [PcMax,Diluted,Pri,Sec] = DilutionMaxPc(r1,v1,cov1,r2,v2,cov2,HBR,params)
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
%    params  - Run parameters for function PcDilution.m [optional]
%               -- params.PriSecScaling: Scale the combined covariance or
%                  each individually (default: 'both'; use 'separate' to 
%                  scale each cov individually and then report the highest 
%                  Pmax, if diluted. 'primary' or 'secondary' can also be 
%                  used to scale only one covariance)
%               -- params.RefineCA: Adjusts the primary 
%                  and secondary states to an actual TCA since
%                  the millisecond precision in TCA inputs often doesn't 
%                  represent the actual TCA of a conjunction.
%                  (default: true)
% Outputs:
%    PcMax   - Maximum Pc value from the combined primary and secondary
%              covariance scaling Pc-dilution analysis
%    Diluted - Integer indicating conjunction is in Pc dilution region for
%              primary covariance scaling or secondary covariance scaling:
%               Diluted =  0 => No dilution for either case
%               Diluted =  1 => Secondary dilution
%               Diluted = 10 => Primary dilution
%               Diluted = 11 => Primary and secondary (or combined) dilution
%    Pri     - Structure holding the primary   Pc-dilution analysis results
%    Sec     - Structure holding the secondary Pc-dilution analysis results
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
%   PcDilution.m
% Subfunctions: None
% MAT-files required: None
%
% See also: none
%
% Last revision: 2025-OCT-15
%
% ----------------- BEGIN CODE -----------------

    % Set up default parameters
    
    if nargin < 8
        params = [];
    end

    % Scale the combined covariance matrix
    if ~isfield(params,'PriSecScaling') || isempty(params.PriSecScaling)
        params.PriSecScaling = 'both';
    end

    % TCA offset correction
    if ~isfield(params,'RefineCA') || isempty(params.RefineCA)
        params.RefineCA = true;
    end

    if strcmpi(params.PriSecScaling,'separate')
        % Analyze Pc-dilution for scaling primary's covariance
        params.PriSecScaling = 'Primary';
        [Pri.PcUnscaled,Pri.Diluted,Pri.PcMax,Pri.SfMax, ...
            Pri.Pc,Pri.Sf,Pri.conv,Pri.iter] = ...
            PcDilution(r1,v1,cov1,r2,v2,cov2,HBR,params);

        % Analyze Pc-dilution for scaling secondary's covariance
        params.PriSecScaling = 'Secondary';
        [Sec.PcUnscaled,Sec.Diluted,Sec.PcMax,Sec.SfMax, ...
            Sec.Pc,Sec.Sf,Sec.conv,Sec.iter] = ...
            PcDilution(r1,v1,cov1,r2,v2,cov2,HBR,params);

        % Combine the results
        DilutedPS = [Pri.Diluted Sec.Diluted];
        PcMaxPS = [Pri.PcMax Sec.PcMax];
        Diluted = 10*Pri.Diluted + Sec.Diluted;
        if Diluted > 0
            PcMax = max(PcMaxPS(DilutedPS));
        else
            PcMax = Pri.PcUnscaled;
        end
    else
        [PS.PcUnscaled,Diluted,PS.PcMax,PS.SfMax, ...
            PS.Pc,PS.Sf,PS.conv,PS.iter] = ...
            PcDilution(r1,v1,cov1,r2,v2,cov2,HBR,params);
        Pri = PS;
        Sec = PS;
        if strcmpi(params.PriSecScaling,'both') && Diluted
            Diluted = 11;
        elseif contains(params.PriSecScaling,'pri','IgnoreCase',true) && Diluted
            Diluted = 10;
        elseif contains(params.PriSecScaling,'sec','IgnoreCase',true) && Diluted
            Diluted = 1;
        elseif ~Diluted
            Diluted = 0;
        end
        if Diluted>0
            PcMax = PS.PcMax;
        else
            PcMax = PS.PcUnscaled;
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
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D.Hall         | 2019-SEP-06 | Initial Development
% S.Es haghi     | 2025-AUG-22 | Adding the option scale the combined
%                                covariance matrix by default. 
% S.Es haghi     | 2025-AUG-25 | Fix default options within params
%                                structure. Fix minor output issues
% S.Es haghi     | 2025-OCT-15 | Set default TCA correction option to true
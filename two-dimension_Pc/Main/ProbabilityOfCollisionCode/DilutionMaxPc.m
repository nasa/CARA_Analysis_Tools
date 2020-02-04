function [PcMax,Diluted,Pri,Sec] = DilutionMaxPc(r1,v1,cov1,r2,v2,cov2,HBR,params)
% DilutionMaxPc - Evaluates collision probability (Pc) dilution and
%                 the associated maximum Pc value resulting from scaling a
%                 conjunction's 3x3 secondary or 3x3 primary covariance.
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
%
% Outputs:
%    PcMax   - Maximum Pc value from the combined primary and secondary
%              covariance scaling Pc-dilution analysis
%    Diluted - Integer indicating conjunction is in Pc dilution region for
%              primary covariance scaling or secondary covariance scaling:
%               Diluted =  0 => No dilution for either case
%               Diluted =  1 => Secondary dilution but no primary dilution
%               Diluted = 10 => Primary dilution but no secondary dilution
%               Diluted = 11 => Primary dilution and secondary dilution
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
% Last revision: 2019-SEP-06
%
% ----------------- BEGIN CODE -----------------

    % Set up default parameters
    
    if nargin < 8
        params = [];
    end

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
%

function [out] = PcCum_PcRem_BootStrap(TActual,PcActual,NpriActual, ...
                                       PcCut,TMission,PcMis_PcRem, ...
                                       NBoot,alpha,ndxBoot)
% PcCum_PcRem_BootStrap - Calculate and estimate the mission cumulative 
% collision probability, and risk remediation probability.
%
% Syntax: [out] = PcCum_PcRem_BootStrap(TActual,PcActual,NpriActual, ...
%                                       PcCut,TMission,PcMis_PcRem, ...
%                                       NBoot,alpha,ndxBoot)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculate and estimate the mission cumulative collision 
% probability, and risk remediation probability, if the PcMis is specified. 
% Otherwise, PcRem needs to be specified. Calculations based on values from
% actual (surrogate) missions and corresponding ratios.
%
% =========================================================================
%
% Input:
%
%   TActual        - Duration of actual (surrogate) missions                              
%
%   PcActual       - Pc values of actual mission conjunctions         [1xN]
%
%   NpriActual     - Number of actual primary missions
%
%   PcCut          - Minimum Pc cutoff 
%
%   TMission       - Duration of mission
%
%   PcMis_PcRem    - 2-element array of mission Pcum and Pcrem        [1x2]
%                    values.May also be a scalar of just mission
%                    Pcum, and Pcrem estimation will be skipped.
%
%   NBoot          - Optional - Number of bootstrap particles. 
%                    (default = 120)
% 
%   alpha          - Optional - Percentile for [alpha, 1-alpha] confidence 
%                    bounds. 
%                    (default = 0.05)
% 
%   ndxBoot        - Optional - Bootstrap realization of Pc   [NPc x NBoot]
%                    values for actual data set time span. 
%                    An array of actual mission Pc indices 
%                    for each bootstrap particle. 
%                    (default = automatically generated)
%
% =========================================================================
%
% Output:
%
% out              - structure containing interpolated median and 
%                    +/- 1-sigma values for Log10 of Prem, achieved 
%                    number of maneuvers, and achieved Pcum
%
% =========================================================================
%
% Dependencies:
%
%   None
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

% Initializations and defaults

Nargin = nargin;

if Nargin < 9; ndxBoot = []; end

if Nargin < 8; alpha = []; end
if isempty(alpha); alpha = 0.05; end

if Nargin < 7; NBoot = []; end
if isempty(NBoot); NBoot = 120; end

% Get the estimate mode from input Pc values
if numel(PcMis_PcRem) == 1
    PcMis_PcRem = [PcMis_PcRem NaN];
end

PcCumMission = PcMis_PcRem(1);
PcRemMission = PcMis_PcRem(2);

PcRemEstimation = isnan(PcRemMission); % Estimate Prem if input is NaN

% Ratios:
%  Trat = mission-duration /actual-duration  and
%  Prat = mission-primaries/actual-primaries

Trat = TMission/TActual;
Prat = 1/NpriActual;

% Sorted Pc array

Pc = PcActual(PcActual >= PcCut);
Pc = sort(Pc);

% Number of events with Pc > PcCut in actual data set time span

NPc = numel(Pc);

% Estimate actual-duration cumulative Pc

PcActu = Pc;
LPsActu = log(1-PcActu);    
LPsCumActu = cumsum(Prat*LPsActu);
out.LPcCumActualInterp = log10(1-exp(LPsCumActu(end)));

% Estimate mission-duration cumulative Pc

LPsCumActu = cumsum(Prat*Trat*LPsActu);
LPcCumActu = log10(1-exp(LPsCumActu));    
out.LPcCumMissionInterp = LPcCumActu(end);

% Calculate the y = Pc vs x = PcCum curve to use later for interpolation
LPcActu = log10(PcActu);
ndx = ~isinf(LPcCumActu);
xxx = LPcCumActu(ndx);
yyy = LPcActu(ndx);
[uxx,nxx] = unique(xxx);
uyy = yyy(nxx);

% Estimate remediation Pc value

if PcRemEstimation

    % Estimate interpolated PcRem given PcCumMission
    LPcCumMission = log10(PcCumMission);
    LPcRemActu = interp1(uxx,uyy, ...
                         LPcCumMission,'linear','extrap');
    out.LPcRemInterp = min(LPcCumMission,LPcRemActu);
    
else
    
    % Estimate interpolated PcCum given PcRemMission
    LPcRemMission = log10(PcRemMission);    
    LPcCumActu = interp1(uyy,uxx, ...
                         LPcRemMission,'linear','extrap');
    out.LPcCumInterp = LPcCumActu;
    
end

% Calculate bootstrap sequence of PcumActual, PcumMission, 
% and PremMission estimates

out.LPcCumActualBoot  = NaN(NBoot,1);
out.LPcCumMissionBoot = NaN(NBoot,1);

out.ndxBoot = NaN(NBoot,NPc);

if PcRemEstimation
    out.LPcRemBoot = NaN(NBoot,1);
else
    out.LPcCumBoot = NaN(NBoot,1);
end

for nb=1:NBoot
    
    % Bootstrap realization of Pc values for actual data set time span
    
    if isempty(ndxBoot)
        ndx = randi(NPc,[NPc 1]);
        % ndx = [randi(NPc-1,[NPc-1 1]); NPc];
    else
        ndx = ndxBoot(nb,:)';
    end
    out.ndxBoot(nb,:) = ndx';
    
    PcBoot = sort(Pc(ndx));
    
    % Estimate actual-duration cumulative Pc
    
    LPsBoot = log(1-PcBoot);    
    LPsCumBoot = cumsum(Prat*LPsBoot);
    out.LPcCumActualBoot(nb) = log10(1-exp(LPsCumBoot(end)));
   
    % Estimate mission-duration cumulative Pc
    
    LPsCumBoot = cumsum(Prat*Trat*LPsBoot);
    
    LPcCumBoot = log10(1-exp(LPsCumBoot));    
    out.LPcCumMissionBoot(nb) = LPcCumBoot(end);
    
    LPcBoot = log10(PcBoot);
    
    % Estimate bootstrap PcRem or PcCum value
    
    if PcRemEstimation
        
        % Estimate interpolated PcRem given PcCumMission
        LPcRBoot = interp_LogPrem(LPcCumBoot,LPcBoot, ...
                             LPcCumMission);
        out.LPcRemBoot(nb) = min(LPcCumMission,LPcRBoot);

    else

        % Estimate interpolated PcCum given PcRemMission
        
        ndx = LPcBoot < LPcRemMission;
        if all(ndx)
            % PcRem > all PcBoot values
            n2 = NPc;
            idx = find(LPcBoot < LPcBoot(n2));
            n1 = idx(end);
        elseif ~any(ndx)
            % PcRem < all PcBoot vaues
            n1 = 1;
            idx = find(LPcBoot > LPcBoot(n1));
            n2 = idx(1);
        else
            ndx = find(ndx);
            n1 = ndx(end);
            idx = find(LPcBoot > LPcBoot(n1));
            n2 = idx(1);
        end
        
        ndx = [n1 n2];
        LPcCBoot = interp1(LPcBoot(ndx),LPcCumBoot(ndx), ...
                           LPcRemMission,'linear','extrap');
        out.LPcCumBoot(nb) = min(0,LPcCBoot);

    end

end

% Calculate bounding indices for bootstrap confidence limits based on input
% value of alpha
alhlf = alpha/2;
n1 = max(floor(alhlf*NBoot),1);
n2 = min(ceil(NBoot-alhlf*NBoot),NBoot);

% Calculate Pc values and 95% bootstrap ranges

out.LPcCumActualBoot = sort(out.LPcCumActualBoot);
out.LPcCumActualMedian = median_sorted(out.LPcCumActualBoot);
out.LPcCumActualRange = [out.LPcCumActualBoot(n1) out.LPcCumActualBoot(n2)];

out.LPcCumMissionBoot = sort(out.LPcCumMissionBoot);
out.LPcCumMissionMedian = median_sorted(out.LPcCumMissionBoot);
out.LPcCumMissionRange = [out.LPcCumMissionBoot(n1) out.LPcCumMissionBoot(n2)];

if PcRemEstimation
    % Estimate interpolated PcRem given PcCumMission
    out.LPcRemBoot = sort(out.LPcRemBoot);
    out.LPcRemMedian = median_sorted(out.LPcRemBoot);
    out.LPcRemRange = [out.LPcRemBoot(n1) out.LPcRemBoot(n2)];
else
    out.LPcCumBoot = sort(out.LPcCumBoot);
    out.LPcCumMedian = median_sorted(out.LPcCumBoot);
    out.LPcCumRange = [out.LPcCumBoot(n1) out.LPcCumBoot(n2)];
end

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
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
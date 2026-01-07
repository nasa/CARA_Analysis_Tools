function [out] = Prem_Nmvr_Bootstrap( ...
    Npriobs,Nobs,Tobs,Pcobs,Tmis,SecCatGrowth,TransMan,Reduction, ...
    Nboot,alpha,lnPcunc,augevent,PremSearch,rndboot)
% Prem_Nmvr_Bootstrap - Calculate the expected number of events in the 
% mission, accounting for various factors
%
% Syntax: [out] = Prem_Nmvr_Bootstrap( ...
%    Npriobs,Nobs,Tobs,Pcobs,Tmis,SecCatGrowth,TransMan,Reduction, ...
%    Nboot,alpha,lnPcunc,augevent,PremSearch,rndboot)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Calculate the expected number of events in the mission, 
% accounting for the mission duration change, the number of primaries 
% combined to create the observed data, and the specified expected growth 
% in the secondary catalog population
%
% =========================================================================
%
% Input:
%
%   Npriobs      - Number of primary surrogates
% 
%   Nobs         - Number of observed conjunctions
% 
%   Tobs         - Number of years spanned by observed set of 
%                  conjunctions
% 
%   Pcobs        - Pc of observed conjunctions                   [Nobs x 1]
% 
%   Tmis         - Mission duration
% 
%   SecCatGrowth - Secondary catalog growth parameter
% 
%   TransMan     - bool indicating Translational (true) or 
%                  rotational (false) RMMs
% 
%   Reduction    - RMM Pc reduction parameter
%
%   Nboot        - Number of bootstrap particles
% 
%   alpha        - parameter to generate confidence intervals. [Nalpha x 1]
%                  May be an array of values to consider
%                  multiple intervals
% 
%   lnPcunc      - Natural log of Pc uncertainty parameter
% 
%   augevent     - Bool indicating whether to consider 
%                  survivalprobability if mission is 
%                  augmented with aworst-case scenario 
%                  conjunction
% 
%   PremSearch   - Array of Prem search values
% 
%   rndboot      - Bootstrap resampling indices. If passed      [Nboot x 1]
%                  in as empty, will be automatically 
%                  generated
%
% =========================================================================
%
% Output:
%
%   out          - structure containing bootstrap indices and
%                  median/confidence interval values for Nmvr and Pcum
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
Nmisexp = Nobs * (Tmis/Tobs) * (1/Npriobs) * SecCatGrowth;
Nmislo  = floor(Nmisexp);
Nmishi  = ceil(Nmisexp);
Nmisrem = rem(Nmisexp,1);

% Generate bootstrap resampling index and normal deviate sequences

if isempty(rndboot)
    
    % Initialize bootstrap index values and variance values to empty cells

    Nndxboot = NaN(Nboot,1);
    rndboot.ndxboot = cell(Nboot,1);
    rndboot.ndvboot = cell(Nboot,1);
    
    % Generate bootstrap index/variance values if there are any observed Pc
    % events
    
    if Nobs > 0

        for nb=1:Nboot
            
            % Calculate number of Pc values for this bootstrap iteration

            if rand < Nmisrem
                Nndxboot(nb) = Nmishi;
            else
                Nndxboot(nb) = Nmislo;
            end
            
            % Sample normal variates
            rndboot.ndvboot{nb} = randn([Nndxboot(nb),1]);

            % Sample indices

            % Unconstrained combinations found by sampling with replacement
            % (most conservative estimation of sampling variations)
            rndboot.ndxboot{nb} = randi(Nobs,[Nndxboot(nb),1]);

        end
        
    end
    
end

out.rndboot = rndboot;

Nalpha = numel(alpha);
alhlf = alpha/2;

nblo = floor(alhlf*Nboot);
nblo(nblo < 1) = 1;

nbhi = ceil(Nboot-alhlf*Nboot);
nbhi(nbhi > Nboot) = Nboot;

% Allocate row arrays for the Pcum and Nmvr results

Nsearch = numel(PremSearch);

out.Pcummn = NaN(1,Nsearch);
out.Nmvrmn = NaN(1,Nsearch);

out.Pcummd = NaN(1,Nsearch);
out.Nmvrmd = NaN(1,Nsearch);

out.Pcumlo = NaN(Nalpha,Nsearch);
out.Pcumhi = NaN(Nalpha,Nsearch);

out.Nmvrlo = NaN(Nalpha,Nsearch);
out.Nmvrhi = NaN(Nalpha,Nsearch);

% Loop over Prem search values

for ns=1:Nsearch
    
    % Initialize Pcum and Nmvr arrays
    
    Pcumrem = zeros(1,Nboot);
    Nmvrrem = zeros(1,Nboot);
    
    % Perform remdiation analysis for all bootstrap realizations
    
    for nb=1:Nboot
        
        % Calculate the mission survival probability
        
        if isempty(rndboot.ndxboot{nb})
            
            % For no mission Pc events, set mission survival prob to one
            
            Psrv = 1;
            
        else
        
            % Extract the bootstrap sequence of mission Pc values
        
            Pcmis = Pcobs(rndboot.ndxboot{nb});
        
            % Introduce logarithmic uncertainty to mission Pc values

            Pcmis = exp(log(Pcmis)+lnPcunc*rndboot.ndvboot{nb});

            % Remediate any mission Pc values exceeding the Prem threshold

            Pcrem = Pcmis; % Initialize Pcrem
            imvr = Pcmis > PremSearch(ns); % Set requiring RMM

            if any(imvr)
                % Count number of remediations
                Nmvrrem(nb) = sum(imvr);
                % Reduce Pc for remediated events
                if TransMan
                    Pcrem(imvr) = Reduction*PremSearch(ns);
                else
                    Pcrem(imvr) = Reduction*Pcmis(imvr);
                end
            end
            
            Psrv = prod(1-Pcrem);

        end
        
        % Calculate the survival probability if the mission is augmented
        % with one worst-case unremediated event with Pc = Prem
        
        if augevent
            Psrv = Psrv*(1-PremSearch(ns));
        end
        
        % Calculate the cumulative Pc
        
        Pcumrem(nb) = 1-Psrv;
        
    end
    
    % Mean values
    
    out.Pcummn(ns) = mean(Pcumrem);
    out.Nmvrmn(ns) = mean(Nmvrrem);

    % Sort remediation Pc values and find the median and 95% range

    Pcumsrt = sort(Pcumrem);
    out.Pcummd(ns) = median_sorted(Pcumsrt);
    for na=1:Nalpha
        out.Pcumlo(na,ns) = Pcumsrt(nblo(na));
        out.Pcumhi(na,ns) = Pcumsrt(nbhi(na));
    end
    
    % Sort number of manuevers and find the median and 95% range
    
    Nmvrsrt = sort(Nmvrrem);
    out.Nmvrmd(ns) = median_sorted(Nmvrsrt);
    for na=1:Nalpha
        out.Nmvrlo(na,ns) = Nmvrsrt(nblo(na));
        out.Nmvrhi(na,ns) = Nmvrsrt(nbhi(na));
    end

    % Make the output Nmvrmd value be the the max of the median and mean
    
    out.Nmvrmd(ns) = max(out.Nmvrmd(ns),out.Nmvrmn(ns));
    
end

return;
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
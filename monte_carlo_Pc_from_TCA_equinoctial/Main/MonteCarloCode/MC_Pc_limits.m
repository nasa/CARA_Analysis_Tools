function [Pcmd,Pclo,Pchi] = MC_Pc_limits(Nc,Nsample,Nsigma)
%
% Calculate a counting probability and Nsigma confidence interval, using
% Matlab's "binofit" function for estimating max. likelihood probabilities
% for Binomial distribution.
%
% =========================================================================
%
% INPUT:
%
%   Nc          = Number of MC counts                  [NxM]
%   Nsample     = Number of MC samples                 [1x1]
%   Nsigma      = Nsigma confidence level(s)           [1x1]
%                 (Negatives indicate -alpha values.)
%
% OUTPUT:
%
%   Pcmd        = Estimated counting probability       [NxM]
%   Pclo        = Lower confidence limit(s)            [NxM]
%   Pchi        = Upper confidence limit(s)            [NxM]
%   
% =========================================================================
%
% REFERENCES:
%
%   See the Wiki page on "Poisson Distribution."
%
% =========================================================================

if (Nsigma < 0)
    alpha = -Nsigma;
else
    alpha = 1-chi2cdf(Nsigma^2,1);
    if any(alpha == 1)
        error('Specified Nsigma value too large. Try Nsigma < 8.3.');
    end
end

if (alpha >= 1) || (alpha < 1e-300)
    error('Alpha value must be between 1e-300 and 1.');
end

% Count the number of elements

S = size(Nc);
N = prod(S);

% Allocate output arrays

Pclo = zeros(S);
Pcmd = Pclo;
Pchi = Pclo;

% Calculate probabilities and confidence intervals

for n=1:N
    
    % Estimate max.like. probability and confidence interval assuming that
    % the MC counts are binomially distributed
    
    [Pcmd(n),Pcci] = binofit(Nc(n),Nsample,alpha); 

    Pclo(n) = Pcci(1);
    Pchi(n) = Pcci(2);

end

return;
end
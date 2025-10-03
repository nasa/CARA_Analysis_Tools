function [zvalue,pvalue] = binomial_prop_test(p1,n1,p2,n2)
% binomial_prop_test - Binomial proportion difference test
%
% Syntax: [zvalue,pvalue] = binomial_prop_test(p1,n1,p2,n2);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
% Implements the binomial proportion difference algorithm to test the null
% hypothesis that two binomial probabilities, p1(actual) and p2(actual),
% are equal. The function uses two estimators for the probilities, p1 and p2.
% One or both can be the result of a Monte Carlo experiment, with the p1
% probability generated with n1 samples, and p2 generated with n2 samples.
% If one (but not both) probability has a sample size of Inf, then
% the function compares that analytically-computed probability against
% the other, Monte Carlo sampled probability.
%
% This function can be used to compare whether two separate Monte Carlo
% runs to compute the probability of collision (Pc) provide similar
% outputs. Alternatively, it can be used to compare an analytically
% computed Pc against a Monte Carlo Pc.
%
% For Pc calculations, the following p-value cutoffs are suggested:
%   1e-3 < pval         = Pc values are similar
%   1e-6 < pval <= 1e-3 = Pc values are somewhat similar
%          pval <= 1e-6 = Pc values are likely dissimilar
%
% References:
%
% NIST Website:
%  https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/diffprop.htm
%
% Original article:
%  Agresti and Caffo (2000), "Simple and Effective Confidence Intervals for
%  Proportions and Differences of Proportions Result From Adding Two
%  Successes and Two Failures", The American Statistician, Vol. 54, No. 4,
%  pp. 280-288.
%
% =========================================================================
%
% Input:
%
%    p1 - Probability of the first sample
%
%    n1 - Number of trials of the first sample (can be Inf, but not if n2
%         is also Inf)
%
%    p2 - Probability of the second sample
%
%    n2 - Number of trials of the second sample (can be Inf, but not if n1
%         is also Inf)
%
% =========================================================================
%
% Output:
%
%    zvalue - Agresti and Caffo z-value
%
%    pvalue - p-value of the comparison
%
% =========================================================================
%
% Initial version: Aug 2023;  Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------

if (n1 <= 0) || (n2 <= 0)
    error('Both n1 and n2 must be positive');
end

if isinf(n1) && isinf(n2)
    error('Both n1 and n2 cannot be infinite');
end

% Calculate Agresti and Caffo (2000) p1-tilde and p2-tilde in a way that
% can tolerate infinite n1 or n2 values
ps1 = (p1+1/n1)/(1+2/n1);
ps2 = (p2+1/n2)/(1+2/n2);

% Calculate Agresti and Caffo (2000) z-value
zvalue = (ps1-ps2)/sqrt(ps1*(1-ps1)/(n1+2)+ps2*(1-ps2)/(n2+2));

% Calculate p-value
[~,pvalue] = ztest(zvalue,0,1);

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
% D. Hall        | 2023-08-22  | Initial Development
% L. Baars       | 2025-09-03  | Updated for public release of code.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

function [Nsamp] = EstimateRequiredSamples(Pc,acc,conf,maxiter)
% EstimateRequiredSamples - Estimates the number of required Monte-Carlo
%                           sample trials needed to receive the accuracy
%                           and confidence levels needed for the Pc passed
%                           in.
%
% Syntax: [Nsamp] = EstimateRequiredSamples(Pc,acc,conf,maxiter);
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
%   Estimate the number of samples required for a Monte Carlo (MC)
%   simulation to achieve a Pc-estimation accuracy of 0 < acc < 1, to a
%   confidence level of 0 < conf < 1.
%
% =========================================================================
%
% INPUT:
%
%   Pc      = Pc to use for estimation process (typically a 2DPc value)
%   acc     = Desired accuracy for the MC Pc estimate (typically 0.1 or
%             0.3)
%   conf    = Confidence level for the accuracy estimate (typically 0.95)
%   maxiter = Maximum number of iterations for binofit analysis [optional]
%              (maxiter <= 0) => Perform no binofit iterations
%              (maxiter >  0) => Perform at most maxiter binofit iterations
%              [default maxiter = 100; typical number of iterations 3 to 10]
%
% =========================================================================
%
% OUTPUT:
%
%   Nsamp   = Recommended number of MC samples to achieve an MC estimation
%             of a probability value 'Pc' with an accuracy of 'acc' at a
%             confidence level of 'conf.'
%
% =========================================================================
%
% Initial version: Mar 2023; Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

% Check for sufficient input

Nargin = nargin;

if Nargin < 3
    error('Insufficient input');
end

% Return infinity for non-positive Pc

if Pc <= 0
    Nsamp = Inf;
    return;
end

% Check for sensible input

if (acc < 0) || (acc > 1)
    error('The input desired accuracy must be 0 < acc < 1, typically 0.1 or 0.3');
end

if (conf < 0) || (conf > 1)
    error('The input desired confidence must be 0 < confidence < 1, typically 0.95 or 0.99');
end

% Set default maximum iterations

if (Nargin < 4) || isempty(maxiter) || (maxiter == Inf)
    maxiter = 100;
end

% Calculate the number of sigmas correspond to the confidence interval (for
% one degree of freedom)

Nsig = sqrt(chi2inv(conf,1));

% Calculate the approximate number of samples to achieve the accuracy at
% the 1-sigma confidence level

Nc1sig = 1/acc^2;
Ns1sig = Nc1sig/Pc;

% Adjust to calculate at the specified confidence interval

Ns0 = max(1,ceil(Nsig^2 * Ns1sig));

Ls = floor(log10(Ns0));
if (Ls > 1)
    Ls = Ls-1;
end

Ms = 10^Ls;
Ns1 = ceil(Ns0/Ms);
Ns = Ns1*Ms;

% Adjust to make binofit return both one-sided accuracies to within the
% specified accuracy at the specified confidence

if (maxiter > 0)
    
    % Calculate the alpha value from the confidence

    alpha = (1-conf);

    % Make the iteration increment for Nsamp approximately 3%
    
    del = max(1,ceil(0.03*Ns1));
    
    % Iterate until binofit indicates the desired accuracy has been
    % achieved, or the maximum iteration number has been exceeded.
    
    niter = 0;
    iterating = true;

    while iterating
        % Calculate binofit confidence interval, Pclo < Pc < Pchi
        [~,Pclohi] = binofit(round(Pc*Ns),Ns,alpha);
        % disp(['Ns = ' num2str(Ns) ' Acc = ' num2str(max(abs(Pclohi-Pc))/Pc)]);
        niter = niter+1;
        % Ensure the estimated accuracy is less than the desired accuracy
        % on both the low and high sides of the interval
        if max(abs(Pclohi-Pc)) > acc*Pc
            % The desired accuracy has not been achieved, so increment
            % the sample number
            Ns = round(Ns/Ms+del)*Ms;
            % Discontinue iterations if max value exceeded
            if (niter >= maxiter)
                iterating = false;
            end
        else
            % The desired accuracy has been achieved, so stop iterating
            iterating = false;
        end
    end
    
end

% Return final estimate for number of samples

Nsamp = Ns;

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
% L. Baars       | 2023-Mar-13 | Initial Development
% L. Baars       | 2025-Aug-06 | Minor documentation updates necessary for
%                                public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
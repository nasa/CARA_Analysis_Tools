function [Nsamp] = EstimateRequiredSamples(Pc,acc,conf,maxiter)

%
% Estimate the number of samples required for a Monte Carlo (MC) simulation
% to achieve a Pc-estimation accuracy of 0 < acc < 1, to a confidence level
% of 0 < conf < 1.
%
% INPUT:
%
% Pc      = Pc to use for estimation process (typically a 2DPc value)
% acc     = Desired accuracy for the MC Pc estimate (typically 0.1 or 0.3)
% conf    = Confidence level for the accuracy estimate (typically 0.95)
% maxiter = Maximum number of iterations for binofit analysis [optional]
%            (maxiter <= 0) => Perform no binofit iterations
%            (maxiter >  0) => Perform at most maxiter binofit iterations
%            [default maxiter = 100; typical number of iterations 3 to 10]
%
% OUTPUT:
%
% Nsamp   = Recommended number of MC samples to achieve an MC estimation
%           of a probability value 'Pc' with an accuracy of 'acc' at a
%           confidence level of 'conf.'
%

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

% Check for valid accuracy and confidence inputs

if (acc <= 0) || (acc >= 1)
    error('The input desired accuracy must be 0 <= acc <= 1, typically 0.1 or 0.3');
end

if (conf <= 0) || (conf >= 1)
    error('The input desired confidence must be 0 <= confidence <= 1, typically 0.95 or 0.99');
end

% Set default maximum iterations

if (Nargin < 4) || isempty(maxiter) || isnan(maxiter)
    maxiter = 100;
elseif (maxiter < 0) || isinf(maxiter)
    error('Invalid maxiter parameter, typically 100');
end

% Calculate the number of sigmas correspond to the confidence interval
% (for a chi-squared distribution with one degree of freedom)

use_statistics_toolbox = license('test','statistics_toolbox');

if use_statistics_toolbox
    % Calculate Nsig using chi2inv
    Nsig = sqrt(chi2inv(conf,1));
else
    % Use zeroth iteration estimate without stats toolbox
    maxiter = 0;
    % Approximate Nsig using table interpolations
    if conf < 1e-2
        % For very small confidences < 0.01
        % (likely never used, but here for completeness)
        Nsig = 1.2533*conf;
    elseif conf <= 9.99e-01
        % For intermediate range confidences: 0.01 to 0.999
        cn = [ ...
            1.0000e-02  1.2533e-02; ...
            1.0000e-01  1.2566e-01; ...
            2.0000e-01  2.5335e-01; ...
            3.0000e-01  3.8532e-01; ...
            4.0000e-01  5.2440e-01; ...
            5.0000e-01  6.7449e-01; ...
            6.0000e-01  8.4162e-01; ...
            6.8269e-01  1.0000e+00; ...
            7.0000e-01  1.0364e+00; ...
            8.0000e-01  1.2816e+00; ...
            8.5000e-01  1.4395e+00; ...
            9.0000e-01  1.6449e+00; ...
            9.3000e-01  1.8119e+00; ...
            9.5000e-01  1.9600e+00; ...
            9.6000e-01  2.0537e+00; ...
            9.7000e-01  2.1701e+00; ...
            9.8000e-01  2.3263e+00; ...
            9.9000e-01  2.5758e+00; ...
            9.9500e-01  2.8070e+00; ...
            9.9900e-01  3.2905e+00];
        Nsig = interp1(cn(:,1),cn(:,2),conf);
    elseif 1-conf >= 1e-11
        % For confidences within 1e-11 to 1e-3 of one
        % (likely never used, but here for completeness)
        c = log10(1-conf);
        cn = [ ...
          -3.0000e+00   3.2905e+00; ...
          -4.0000e+00   3.8906e+00; ...
          -5.0000e+00   4.4172e+00; ...
          -6.0000e+00   4.8916e+00; ...
          -7.0000e+00   5.3267e+00; ...
          -8.0000e+00   5.7307e+00; ...
          -9.0000e+00   6.1094e+00; ...
          -1.0000e+01   6.4670e+00; ...
          -1.1000e+01   6.8065e+00];
        Nsig = interp1(cn(:,1),cn(:,2),c);
    else
        error('Confidence out of tabulated range');
    end
end

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
        if use_statistics_toolbox
            [~,Pclohi] = binofit(round(Pc*Ns),Ns,alpha);
        else
            Pclohi = sqrt(round(Pc*Ns))/Ns;
        end
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
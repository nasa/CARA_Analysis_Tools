
% determines a single scale factor for the covariance that minimizes one of selected goodness-of-fit test
% statistics for testing for a normalized residual set's confirmity to a 3-DoF chi-squared distribution.
% Two approaches are taken here:  an optimization using the entire residual set, and a resampled approach
% that outputs a PDF of scale factors.  While it could be done from the results of this routine, the
% intent is for the user to take the optimized scale factor and re-run the analysis in
% EvaluateResiduals.m to get a more complete graphical presentation of the performance with that scale
% factor employed.  Furthermore, it is not really recommended to use scale factors in operations, and
% certainly not based on a short analysis as outlined here; rather, the scale factor gives a sense of
% whether the covariance is too large or small, and approximately by how much.
% -------------------------------------------------------------------------------------------------------
% Inputs:
%   Residuals:      n x 3 set of position residuals
%   Covariances:    n x 3 x 3 set of position covariances matching the residuals
%   Options:        (optional) controls a variety of options for the covariance realism analysis:
%                   NumberOfTrials:  number of sample trials to take from dataset for GOF testing (10000 recommended/default)
%                   TrialSampleSize:  number of samples from dataset to take for each trial  (50 recommended/default)
%                   SignificanceLevel:  GOF test p-value (0.02 recommended)
%   TestStat:       (optional) the statistic to use when calculating optimal scale factors (1 - Cramer - von Mises; 3 =
%                   Anderson-Darling) (Default = 1)
%   FigureOffset:   (optional) offsets to use in applying Figure numbers; useful if analyzing several different.
%                   Recommended value = 10; set to 0 to suppress generation
%                   of graphs (default = 0)
% -------------------------------------------------------------------------------------------------------
% Outputs:
%   FullM2ScaleFactor:      scale factor, based on the full set of residuals, that minimizes the selected
%                           test statistic
%                           Also, a graph is produced that gives a PDF of the scale factors observed during the resampling
%                           analysis
% -------------------------------------------------------------------------------------------------------
%   .m file dependencies:
%       FullSpecDistVec.m
%       labelpoints.m
% -------------------------------------------------------------------------------------------------------
%   Developed by M.D. Hejduk in support of NASA/CARA
%   Copyright © 2020 United States Government as represented by the Administrator of the National 
%   Aeronautics and Space Administration.  All Rights Reserved.

function FullM2ScaleFactor=DetermineOptimalScaleFactor(...
        Residuals,Covariances,Options,TestStat,FigureOffset)
    %% Set input values if not defined in inputs
    Nargin = nargin;
    
    % Populate undefined Options
    if Nargin < 3 || isempty(Options)
        Options = struct();
    end
    if ~isfield(Options,'NumberOfTrials') || isempty(Options.NumberOfTrials)
        Options.NumberOfTrials=10000;
    end
    if ~isfield(Options,'TrialSampleSize') || isempty(Options.TrialSampleSize)
        Options.TrialSampleSize=50;
    end
    if ~isfield(Options,'SignificanceLevel') || isempty(Options.SignificanceLevel)
        Options.SignificanceLevel=0.02;
    end
    
    
    if Nargin < 4 || isempty(TestStat)
        TestStat = 1;
    end
    
    % If no figure offset defined, suppress Output
    if Nargin < 5 || isempty(FigureOffset)
        FigureOffset = 0;
    end
    
    
    
    %% Analyze Inputs
    n=size(Residuals,1);
    % creates (squares of) Mahalanobis distances
    M2=zeros(n,1);
    for i=1:1:n
        % square of Mahalanobis distance:  e * inv(C) * eT, in which e is a single residual vector (1 x 3
        % here) and C is the covariance (3 x 3) associated with the residual.  The backslash operator is
        % used instead of the inv() function because it is both more stable and more efficient
        M2(i)=Residuals(i,:)/squeeze(Covariances(:,:,i))*Residuals(i,:)';
    end
    % initial guess of scale factor (set to unity)
    s0=1;
    % anonymous function that calculates the desired test statistic for a particular scale factor
    % setting.  The function takes as inputs the vector of Mahalanobis distances and a flag for which
    % test statistic to evaluate:
    %   1 = Cramer - von Mises statistic
    %   2 = Watson statistic (included because of some legacy code; DO NOT use this statistic here)
    %   3 = Anderson-Darling statistic
    % finds the optimal (square of) scale factor, using the test statistic indicated, against the entire
    % residual set
    fun=@(s)TransformAndCalculate(s,M2,TestStat);
    [SF,fminval,exitflag]=fminsearch(fun,s0);
    % more elaborate preservation of results from minimization search
    FullM2ScaleFactorSet=[sqrt(SF) fminval exitflag];
    % preserves just the scale factor for export
    FullM2ScaleFactor=FullM2ScaleFactorSet(1);
    % preallocation of "square of scale factor" structure.  Typically (and in this implementation) the
    % scaling of the covariance is through pre- and post-multiplication, so covariances are effectively
    % scaled by the square of the scale factor
    ScaleFactorsSquared=NaN(Options.NumberOfTrials,3);
    % resampling approach to calculate a PDF of scale factors.  This is a good way to address the
    % deleterious effects of outliers on the GOF testing, and it establishes a reasonable range of scale
    % factors and thus shows how tightly or broadly they are distributed
    parfor j=1:Options.NumberOfTrials
        % generates the random samples (MATLAB function here does not support vectorization)
        Samples=randsample(M2,Options.TrialSampleSize);
        % initial guess at scale factor; nonlinear minimization routine requires this.  Since the test
        % statistic is of the order y=f(x)^2, the curve should be at least quasi-parabolic; so the initial
        % guess should not matter--the problem should be extremely well-behaved and find the global minimum
        % quickly.  fminsearch is probably overkill for this problem, actually.
        s0=1;
        % anonymous function that calculates the desired test statistic for a particular scale factor
        % setting.  The function takes as inputs the vector of Mahalanobis distances and a flag for which
        % test statistic to evaluate:
        %   1 = Cramer - von Mises statistic
        %   2 = Watson statistic (included because of some legacy code; DO NOT use this statistic here)
        %   3 = Anderson-Darling statistic
        % finds the optimal scale factor, using the test statistic indicated
        fun=@(s)TransformAndCalculate(s,Samples,TestStat);
        [s1,fminval,exitflag]=fminsearch(fun,s0);
        ScaleFactorsSquared(j,:)=[s1 fminval exitflag];
    end
    % square root of resampled scale factors 
    ScaleFactors=sqrt(ScaleFactorsSquared(:,1));
    if FigureOffset>0
        figure(FigureOffset+1);
        % synthesizes PDF for output
        LowMedHigh=prctile(ScaleFactors,[1 50 99]);
        pd=fitdist(ScaleFactors,'kernel','Kernel','epanechnikov');
        X=LowMedHigh(1):0.001:LowMedHigh(3);
        Y=pdf(pd,X);
        % places scale factor for full M2 set onto resampled PDF
        YForSF=interp1(X,Y,FullM2ScaleFactor);
        h=plot(X,Y,FullM2ScaleFactor,YForSF,'*m');
        set(h(2),'linewidth',2);
        grid on
        legend('PDF Estimate','Full Sample Scale Factor','location','south');
        xlabel('Scale Factor')
        set(gca,'yticklabel',[]);
        title('Scale Factor Estimated PDF');
        labelpoints(FullM2ScaleFactor,YForSF,num2str(FullM2ScaleFactor),'NW');
    end
end



% function to calculate the test statistic for each iteration of the scale factor solution routine.
% Rather than run the optimizer with the p-value output (which could fail if the results push outside of
% the range of the interpolatable tables for p-value calculation), it is more robust to work directly
% from the test statistic itself and minimize that.  
% The scale factors are in principle applied to the covariances themselves, but one can apply them to
% Mahalanobis distance histories by dividing them by the square of the scale factors
%
% Inputs:
%   x:  current value of squared scale factor (in current iteration of minimization routine)
%   M2: unscaled Mahalanobis distance history
%   TestStat:  desired test statistic to use (1:  Cramer - von Mises; 3:  Anderson-Darling)
% Outputs:
%   Output:  value of selected test statistic

function Output=TransformAndCalculate(x,M2,TestStat)
    % x square of scale factor
    M2Norm=M2./x;
    % call to goodness-of-fit test routine
    [~,QTrans]=FullSpecDistVec(M2Norm,3,'chi2');
    % value of selected test statistic for this iteration
    Output=QTrans(TestStat);
end


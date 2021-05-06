
% simple driver program for covariance realism analysis
function AnalysisDriver(FileName)
    % as explained elsewhere, resampling will be used to lessen the effect of outliers.  Here the number
    % of trials to execute is specified.  For anything definitive, at least 10,000 trials should be used.
    % The analysis results for the entire sample are always returned
    Options.NumberOfTrials=5000;
    % the performance of the goodness-of-fit (GOF) tests are, unfortunately, affected by sample size,
    % despite normalization techniques. A size of around 50 is probably good; sizes substantially smaller
    % will bias the results in a favorable direction.
    Options.TrialSampleSize=50;
    % the p-value to use for the GOF test results is also a matter of preference.  0.05 is a "standard"
    % of sorts; for covariance realism applications, going down to perhaps 0.02 would be acceptable; but
    % a lower value than that is unwise.
    Options.SignificanceLevel=0.02;
    % allows graphs to display in their title the propagation state of the data analyzed.  
    PropagationStateText='x Days';
    % loads data to test for realism
%     load('TestInput01');
%     load('TestInput02');
%     load('TestInput03');
%     load('TestInput04');
%     load('TestInput05');
    load(FileName);
    % sets the scale factor to use for a particular evaluation.  For the first evaluation, the scale
    % factor is set to unity so that the unscaled covariance realism performance can be evaluated.
    ScaleFactor=1;
    % keeps generated graphs from overwriting each other.  Approach here is to use a sequential multiple
    % of 10 for each sequenced program run
    FigureOffset=10;
    % initial evaluation of dataset:  will produce graphs evaluating normality for each component errors
    % and chi-squared compliance for the squared Mahalanobis distance distribution
    [NormalityTestFullResiduals,NormalityTestResampledResidualsInterpolatable, ...
    Chi2TestFullM2,Chi2TestResampledM2Interpolatable]=EvaluateResiduals(Residuals,Covariances,ScaleFactor,Options,FigureOffset,PropagationStateText,1);
    % determines which testing approach to be used for the scale factor evaluation.  1 = Cramer - von
    % Mises; 3 = Anderson-Darling (2 is not used in this context)
    TestStat=1;
    FigureOffset=20;
    % evaluates scale factors that could be used to make covariance results better.  These are not
    % intended to be used as an actual operational fix (i.e., to multiply the generated operational
    % covariance by this scale factor, but rather to indicate the general level of covariance mis-sizing
    % and therefore how to proceed with tuning
    FullSetScaleFactor=DetermineOptimalScaleFactor(Residuals,Covariances,Options,TestStat,FigureOffset);
    FigureOffset=30;
    % re-performs the original analysis but with the covariances scaled by the scale factor determined
    % above
    EvaluateResiduals(Residuals,Covariances,FullSetScaleFactor,Options,FigureOffset,PropagationStateText,0);
end
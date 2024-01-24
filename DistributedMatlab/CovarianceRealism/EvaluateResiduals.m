
% Evaluates normalized residuals both for Gaussian behavior by individual component and 2-DoF chi-squared
% behavior as ordered triples of position information.  The normalization is by the covariance associated
% with each residual ordered triple, and is performed within this routine.  Evaluation is through an
% empirical distribution function goodness-of-fit (GOF) technique, either the Cramer - von Mises or
% Anderson-Darling test.
% -------------------------------------------------------------------------------------------------------
% Inputs:
%   Residuals:      n x 3 set of position residuals
%   Covariances:    n x 3 x 3 set of position covariances matching the residuals
%   Scale factor:   (optional) scale factor by which to multiply the covariance.  Because a pre- and
%                   post-multiplication is presumed, covariance is actually multiplied by the square 
%                   of the scale factor (default = 1)
%   Options:        (optional) controls a variety of options for the covariance realism analysis:
%                   NumberOfTrials:  number of sample trials to take from dataset for GOF testing (10000 recommended/default)
%                   TrialSampleSize:  number of samples from dataset to take for each trial  (50 recommended/default)
%                   SignificanceLevel:  GOF test p-value (0.02 recommended/default)
%   FigureOffset:   (optional) offsets to use in applying Figure numbers; useful if analyzing several different.
%                   Recommended value = 10; set to 0 to suppress generation
%                   of graphs (Default = 0)
%   PropagationStateText:  (optional) allows the propagation state to be
%                           displayed in graph titles (default = '')
%   MeansGraph:     (optional, Default = false) flag to determine whether to produce graph of means (needed only for first-pass
% -------------------------------------------------------------------------------------------------------
% Outputs:
%   NormalityTestFullResiduals:     3 x 4 results from testing for Gaussian behavior the entire residual
%                                   set, by component.  Each row is the results from a particular test 
%                                   (Cramer - von Mises, Watson [not used], and Anderson-Darling), and 
%                                   each column is a component (radial, in-track, cross-track, and a
%                                   blank column to keep the array from being 3 x 3 and thus confusing).  
%                                   Array gives p-value:  0.009 means < 0.01; 0.26 means > 0.25 
%   NormalityTestResampledResidualsInterpolatable:  same idea as the previous array, but here resampled
%                                   results are given (percentage of trials that passed, according to 
%                                   specified p-value).  Same definition of array contents 
%   Chi2TestFullM2:                 gives results of the test of full normalized position error vectors (full dataset)
%                                   for conformity to a 3-DoF chi-squared distribution.  Array is 3 x 4, but only 
%                                   column 4 is populated. Array contains p-value information
%   Chi2TestResampledM2Interpolatable:  same as previous array but with results from all of the trials
%                                   from the resampling approach.  Same array definitions, but percent of 
%                                   cases passing (at the specified p-value) is what is reported
% -------------------------------------------------------------------------------------------------------
%   .m file dependencies:
%       FullSpecDistVec.m
% -------------------------------------------------------------------------------------------------------
%   Developed by M.D. Hejduk in support of NASA/CARA
%   Copyright © 2020 United States Government as represented by the Administrator of the National 
%   Aeronautics and Space Administration.  All Rights Reserved.

function [NormalityTestFullResiduals,NormalityTestResampledResidualsInterpolatable, ...
    Chi2TestFullM2,Chi2TestResampledM2Interpolatable]=EvaluateResiduals(...
        Residuals,Covariances,ScaleFactor,Options,FigureOffset,PropagationStateText,MeansGraph)
    %% Set input values if not defined in inputs
    Nargin = nargin;
    
    % Populate Scale Factor
    if Nargin < 3 || isempty(ScaleFactor)
        ScaleFactor = 1;
    end
    
    % Populate undefined Options
    if Nargin < 4 || isempty(Options)
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
    
    % If no figure offset defined, suppress Output
    if Nargin < 5 || isempty(FigureOffset)
        FigureOffset = 0;
    end
    
    if Nargin < 6 || isempty(PropagationStateText)
        PropagationStateText = '';
    end
    
    if Nargin < 7 || isempty(MeansGraph)
        MeansGraph = false;
    end
    
    %% Analyze Inputs
    
    % applies scale factor to covariances (note square of scale factor, since process is actually a pre-
    % and post-multiplication)
    Covariances=Covariances*ScaleFactor^2;
    
    % number of residuals sets and covariances to evaluate (for a fixed propagation time)
    n=size(Residuals,1);
    % ---------------------------------------------------------------------------------------------------
    % tests each components' normalized residuals for normality
    % structure for normalized residuals.  Purpose here is not so much to check for residual Gaussian
    % behavior but whether the residuals normalized by the covariance standard deviations conform to a
    % standardized Gaussian distribution
    NormResid=NaN(n,3);
    % output structure for normality test of full set of each component's normalized residuals
    NormalityTestFullResiduals=zeros(3,4);
    % output structure for normality test of resampled sets of each component's normalized residuals
    NormalityTestResampledResiduals=zeros(3,Options.NumberOfTrials,3);
    % resampled normality results expressed as a CDF
    NormalityTestResampledResidualsInterpolatable=zeros(3,4);
    % executes once for each 
    for i=1:1:3
        NormResid(:,i)=Residuals(:,i)./sqrt(squeeze(Covariances(i,i,:)));
        % calls distribution testing routine, which will invoke both Cramer - von Mises and
        % Anderson-Darling tests for a fully-specified distribution (zero-mean, unity variance normal);
        % here fthe entire dataset is tested
        NormalityTestFullResiduals(:,i)=FullSpecDistVec(NormResid(:,i),[0 1]','normal');
        % preallocates structure for random samples from residuals (to test each sample for normality)
        Samples=zeros(Options.TrialSampleSize,Options.NumberOfTrials);
        % generates the random samples (MATLAB function here does not support vectorization)
        for j=1:1:Options.NumberOfTrials
            Samples(:,j)=randsample(NormResid(:,i),Options.TrialSampleSize);
        end
        % calculates the mean value for each sample set
        TheMeans=mean(Samples);
        % synthesizes PDF for output
        LowMedHigh=prctile(TheMeans,[1 50 99]);
        pd=fitdist(TheMeans','kernel','Kernel','epanechnikov');
        Means{i}(:,1)=LowMedHigh(1):0.001:LowMedHigh(3);
        Means{i}(:,2)=pdf(pd,Means{i}(:,1));
        % calls distribution testing routine, which will invoke both Cramer - von Mises and
        % Anderson-Darling tests for a fully-specified distribution (zero-mean, unity variance normal);
        % here, one sample from the full dataset is tested
        NormalityTestResampledResiduals(:,:,i)=FullSpecDistVec(...
            Samples,repmat([0 1]',1,Options.NumberOfTrials),'normal');
        % creates CDF plots of p-values for each sample group tested
        [f{i,1},x{i,1}]=ecdf(NormalityTestResampledResiduals(1,:,i));
        % creates CDF output of resampled results for Cramer - von Mises test
        NormalityTestResampledResidualsInterpolatable(1,i)=OutputInterpolation(...
            x{i,1},f{i,1},Options.SignificanceLevel);
        % creates CDF output of resampled results for Anderson-Darling test
        [f{i,3},x{i,3}]=ecdf(NormalityTestResampledResiduals(3,:,i));
        NormalityTestResampledResidualsInterpolatable(3,i)=OutputInterpolation(...
            x{i,3},f{i,3},Options.SignificanceLevel);
    end
    %% if FigureOffset greater than 0, generates output plots
    % plotting routine labeling presumes residuals are "radial, in-track, cross-track"
    if FigureOffset>0
        % creates pdf plots of resampled means (if MeanGraph set to 1)
        if MeansGraph
            figure(FigureOffset+1)
            plot(Means{1}(:,1),Means{1}(:,2),Means{2}(:,1),Means{2}(:,2),Means{3}(:,1),Means{3}(:,2));
            grid on
            xlabel('Mean of Normalized Residual')
            legend('Radial','In-Track','Cross-Track','location','northeast');
            title('Estimated PDFs of Resampled Normalized Residuals');
        else
            FigureOffset=FigureOffset-1;
        end
        % creates plots of normalized component residuals' Gaussian distribution compliance
        % Cramer - von Mises test
        figure(FigureOffset+2)
        subplot(1,2,1)
        plot(x{1,1},100-f{1,1}*100,x{2,1},100-f{2,1}*100,x{3,1},100-f{3,1}*100);
        grid on
        legend('Radial','In-Track','Cross-Track','location','northeast');
        xlabel('p-Value');
        ylabel('Cumulative Percentage');
        title('Cramer - von Mises Test')
        set(gca,'xlim',[0.01 0.25]);
        set(gca,'xtick',[0.01 0.05 0.1 0.15 0.2 0.25],'xticklabel',[0.01 0.05 0.1 0.15 0.2 0.25]);
        % Anderson-Darling test
        subplot(1,2,2)
        plot(x{1,3},100-f{1,3}*100,x{2,3},100-f{2,3}*100,x{3,3},100-f{3,3}*100);
        grid on
        legend('Radial','In-Track','Cross-Track','location','northeast');
        xlabel('p-Value');
        ylabel('Cumulative Percentage');
        title('Anderson-Darling Test');
        set(gca,'xlim',[0.01 0.25]);
        set(gca,'xtick',[0.01 0.05 0.1 0.15 0.2 0.25],'xticklabel',[0.01 0.05 0.1 0.15 0.2 0.25]);
        if ~isempty(PropagationStateText)
            sgtitle(['Portion of Samples Passing Normality Test:  ' PropagationStateText ...
                ', Scale Factor = ' num2str(ScaleFactor)]);
        else
            sgtitle(['Portion of Samples Passing Normality Test: Scale Factor = ' num2str(ScaleFactor)]);
        end
    end
    % ---------------------------------------------------------------------------------------------------
    %% Test the residuals' (squares of) Mahalanobis distances for conformity to a 3-DoF chi-square
    % distribution 
    % creates (squares of) Mahalanobis distances
    M2=zeros(n,1);
    for i=1:1:n
        % square of Mahalanobis distance:  e * inv(C) * eT; e is vector of residuals and C covariance
        M2(i)=Residuals(i,:)/squeeze(Covariances(:,:,i))*Residuals(i,:)';
    end
    % calls distribution testing routine, which will invoke both Cramer - von Mises and
    % Anderson-Darling tests for a fully-specified distribution (3-DoF chi-squared distribution)
    % first, the entire distribution is tested
    Chi2TestFullM2(:,4)=FullSpecDistVec(M2,3,'chi2');
    % preallocates structure for random samples 
    Samples=zeros(Options.TrialSampleSize,Options.NumberOfTrials);
    % generates the random samples (MATLAB function here does not support vectorization)
    for j=1:1:Options.NumberOfTrials
        Samples(:,j)=randsample(M2,Options.TrialSampleSize);
    end
    % calls distribution testing routine, which will invoke both Cramer - von Mises and
    % Anderson-Darling tests for a fully-specified distribution (3-DoF chi-squared distribution)
    Chi2TestResampledM2(:,:,4)=FullSpecDistVec(...
        Samples,ones(1,Options.NumberOfTrials)*3,'chi2');
    % creates CDF plots of p-values for Cramer - von Mises test
    [f{4,1},x{4,1}]=ecdf(Chi2TestResampledM2(1,:,4));
    Chi2TestResampledM2Interpolatable(1,4)=OutputInterpolation(x{4,1},f{4,1},Options.SignificanceLevel);
    %% Create CDF plots of p-values for Anderson-Darling test
    [f{4,3},x{4,3}]=ecdf(Chi2TestResampledM2(3,:,4));
    Chi2TestResampledM2Interpolatable(3,4)=OutputInterpolation(x{4,3},f{4,3},Options.SignificanceLevel);
    if FigureOffset>0
        figure(FigureOffset+3)
        plot(x{4,1},100-f{4,1}*100,x{4,3},100-f{4,3}*100);
        grid on
        set(gca,'xlim',[0.01 0.25]);
        set(gca,'xtick',[0.01 0.05 0.1 0.15 0.2 0.25],'xticklabel',[0.01 0.05 0.1 0.15 0.2 0.25]);
        xlabel('p-Value');
        ylabel('Cumulative Percentage');
        legend('Cramer - von Mises Test','Anderson-Darling Test','location','northeast');
        if ~isempty(PropagationStateText)
            title(['Mahalanobis Distance Conformity to 3-DoF \chi^2 Distribution:  ' ...
                PropagationStateText ', Scale Factor = ' num2str(ScaleFactor)]);
        else
            title(['Mahalanobis Distance Conformity to 3-DoF \chi^2 Distribution: Scale Factor = ' num2str(ScaleFactor)]);
        end
    end
end


% MATLAB edcf function can include duplicate x-values (for stair-step behavior), and this breaks
% interpolators.  This routine adds an insignificant sequential small value to each x-value so that all
% x-values will be unique, allowing interpolation routines to work correctly; and then it calls the
% one-dimensional interpolator
function Output=OutputInterpolation(x,f,pValue)
    n=length(x);
    Offset=((1:1:n)*0.0000001)';
    % interpolation routine.  f multiplied by 100 to allow display as regular percentages
    Output=interp1(x+Offset,f*100,pValue,'pchip');
end














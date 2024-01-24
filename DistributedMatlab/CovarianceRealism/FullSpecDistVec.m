
% performs an Empirical Distribution Function (EDF) goodness-of-fit test for a fully-specified
% distribution (the so-called "case 0").  This procedure follows D'Agnostino and Stephens, Goodness of
% Fit Techniques (1986), pp. 97 - 105; users should consult this monograph for more information.  This
% routine is vectorized, allowing multiple sets (which must be of the same size) to be tested as part of
% a single function call
% -------------------------------------------------------------------------------------------------------
% Inputs:
%    x - the data to be tested.  Each column is a separate set to be tested.  All columns must be the
%       same length; if this is not true, successive calls with only one set of data must be used.
%    param - the array of distribution parameters to be used.  The param array is 3 x n, where n is the
%       number of columns in x.  The first row contains the location parameters, the second row the scale
%       parameters, and the third row the shape parameters.  A parameter set must be specified for each
%       column of x
%    DistType - indicates the distribution to be tested.  The different possibilities are:
%       normal - normal distribution (param row 1 is mean, row 2 is standard deviation)
%       exponential - exponential distribution (param row 1 is mean)
%       extremevalue - extreme value distribution (param row 1 is location parameter, 2 is scale param)
%       gamma2p - two-parameter gamma distribution (param row 1 not used, 2 is scale param, 3 is shape)
%       chi2 - chi-square distribution (param row 1 is degrees of freedom)
%       weibull2p - two-parameter Weibull dist (param row 1 not used, 2 is scale param, 3 shape param)
%       weibull3p - three-param Weibull dist (param row 1 location, 2 scale, 3 shape)
% -------------------------------------------------------------------------------------------------------
% Outputs:
%   FullSpecOut is a 3 x n array (n number of columns in x) in which
%       row 1 is the set of p-values for the Cramer - von Mises EDF statistic (W2)
%       row 2 is the set of p-values for the Watson EDF statistic (U2)
%       row 3 is the set of p-values for the Anderson-Darling EDF statistic (A2)
%   For this test, the p-value table runs from 0.001 to 0.25.  Returns of 0.0009 or 0.26 indicate that
%   the test value lies outside of the range of values represented in the table, at the indicated extreme
%   QTrans is same dimension as FullSpecOut but gives the raw test statistic (transformed Q).  This is
%   sometimes helpful in comparing results from different datasets that fall outside of the p-value
%   limits.  Larger is worse.
% -------------------------------------------------------------------------------------------------------
%   .m file dependencies:
%       CalculateEDFVec.m
% -------------------------------------------------------------------------------------------------------
%   Developed by M.D. Hejduk in support of NASA/CARA
%   Copyright © 2020 United States Government as represented by the Administrator of the National 
%   Aeronautics and Space Administration.  All Rights Reserved.

function [FullSpecOut,QTrans]=FullSpecDistVec(x,param,DistType)
    % ensures that input data be sorted
    x=sort(x,1);
    % creates structures of distribution parameters that match the expected different arrays and array
    % sizes exppected by the MATLAB cdf functions
    Param1=repmat(param(1,:),size(x,1),1);
    if size(param,1)>=2
        Param2=repmat(param(2,:),size(x,1),1);
    end
    if size(param,1)==3
        Param3=repmat(param(3,:),size(x,1),1);
    end
    % creates cdf of idealized distribution against which the empirical data in x will be tested; all
    % supported distribution types are shown below.  While it is believed that this software works
    % properly for all of the types shown here, only the normal, exponential, and chi-squared have been
    % rigorously tested
    switch DistType
        case 'normal'
            Z=normcdf(x,Param1,Param2);
        case 'exponential'
            Z=expcdf((x-Param1),Param2);
        case 'extremevalue'
            Z=evcdf(x,Param1,Param2);
        case 'gamma2p'
            Z=gamcdf(x,Param3,Param2);
        case 'chi2'
            Z=chi2cdf(x,Param1);
        case 'weibull2p'
            Z=wblcdf(x,Param2,Param3);
        case 'weibull3p'
            Z=1-exp(-((x-Param1)./Param2).^Param3);
        otherwise
            error('FullSpecDist:  Distribution type not supported');
    end
    % calculates W2, U2, and A2 from the comparison of empirical and ideal CDFs; Q gives the test
    % statistic output 
    Q=CalculateEDFVec(Z);
    % assigns a p-value (significance level) to the W2, U2, and A2 values 
    [FullSpecOut,QTrans]=EvalEDF(Q,size(x,1));
end


% calculates the EDF statistic for the three main EDF test types (Cramer - von Mises, Watson, and
% Anderson-Darling).  Follows the formulae in D'Agostino and Stephens (1986), p. 101
function QStat=CalculateEDFVec(z)
    % number of data points
    n=size(z,1);
    % number of datasets
    m=size(z,2);
    % Cramer - von Mises
    Id=repmat((1:1:n)',1,m);
    I=(2.*Id-1)./(2*n);
    W2=sum((z-I).^2)+1/(12*n);
    % Watson
    U2=W2-n.*(mean(z)-0.5).^2;
    % Anderson-Darling
    Term1=(2.*Id-1).*log(z);
    Term2=(2*n+1-2*Id).*log(1-z);
    A2=-n-(1/n).*sum((Term1+Term2));
    QStat=[W2; U2; A2];
end


% Assigns p-values (significance levels) to the W2, U2, and A2 values calculated by the CalculateEDF
% function.  Q is the test statistic output from CalculateEDF, and n is the number of sample points
% (which should be the same for each results set tested)
function [pValResults,QTrans]=EvalEDF(Q,n)
    % preallocation of variables
    pValResults=zeros(size(Q)); QTrans=pValResults; EDFTable = zeros(3,8);
    % test statistic transformations in order to adjust for sample size.  This turns W2, U2, and A2 into
    % W*, U*, and A*
    QTrans(1,:) = (Q(1,:) - 0.4/n + 0.6/n^2) .* (1 + 1/n);
    QTrans(2,:) = (Q(2,:) - 0.1/n + 0.1/n^2) .* (1 + 0.8/n);
    QTrans(3,:) = Q(3,:);
    % tables giving p-values for test statistic values
    % p-value definitions
    EDFBoundary = [.25 .15 .10 .05 .025 .01 .005 .001];
    % W2 test statistic values corresponding to p-value definitions
    EDFTable(1,:) = [0.209 0.284 0.347 0.461 0.581 0.743 0.869 1.167];
    % U2 test statistic values corresponding to p-value definitions
    EDFTable(2,:) = [0.105 0.131 0.152 0.187 0.222 0.268 0.304 0.385];
    % A2 test statistic values corresponding to p-value definitions
    EDFTable(3,:) = [1.248 1.610 1.933 2.492 3.070 3.880 4.500 6.000];
    % performs interpolation (or boundary indication) on p-value tables to determine p-values for each
    % test statistic
    % one execution for each EDF parameter (W*, U*, A*)
    for i = 1:1:3
        % entries greater than smallest table value; output set to 0.26 as high-end flag
        Index1=QTrans(i,:)<EDFTable(i,1);
        pValResults(i,Index1)=0.26;
        % entries smaller than largest table value; output set to 0.009 as high-end flag
        Index2=QTrans(i,:)>EDFTable(i,end);
        pValResults(i,Index2)=0.0009;
        % interpolation
        pValResults(i,~(Index1 | Index2))=...
            interp1(EDFTable(i,:),EDFBoundary,QTrans(i,~(Index1 | Index2)));
    end
end
    



function ODQualityThresholds=SetODQualityAssessmentThresholds
%
% ODQualityThresholds - Sets Default OD quality assessment thresholds for
%                       use within the "CDM_ODQualityAssessment" routine
%
% Syntax:               ODQualityThresholds=SetODQualityAssessmentThresholds
%
% Inputs: None
%
% Outputs:
%   ODQualityThresholds -   Structure of default threshold values to be
%                           used for Orbit Determination Quality
%                           Assessmenent purposes
%
% Examples/Validation Cases:
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: CDM_ODQualityAssessment.m
%
% July 2019; Last Revision: 24-07-2019
%
% ----------------- BEGIN CODE -----------------
    ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable = ...
                                [0.00   0.25      0    500   36   1   0;
                                 0.00   0.25    500    900   36   1   1;
                                 0.00   0.25    900   2000   24   1   1;
                                 0.25   1.00      0    500   36   1   0;
                                 0.25   1.00    500   1000   24   1   1;
                                 0.25   1.00   1000   2000   18   0   1;
                                 0.00   1.00   2000  10000   12   0   1;
                                 0.00   1.00  10000 100000    8   0   1];
    ODQualityThresholds.BReasonabilityEDR1 =    [0.001 0.1
                                                 0.001 0.2
                                                 0.001 1.0];
    ODQualityThresholds.BReasonabilityEDRGT1 =  [0.001 0.1
                                                 0.001 0.2
                                                 0.001 10];
    ODQualityThresholds.SRPReasonabilityEDR01 = [0.001 0.1
                                                 0.001 0.2
                                                 0.001 1.0];
    ODQualityThresholds.SRPReasonabilityEDRGT1 =[0.001 0.1
                                                 0.001 0.2
                                                 0.001 1.0];
    ODQualityThresholds.ODResults.PercentResidualAcceptance =   80;
    ODQualityThresholds.ODResults.MinWRMS =                     [0 0 0];
    ODQualityThresholds.ODResults.WRMS =                        [1.5 2.0 5.0];
    ODQualityThresholds.ODResults.LUPIRatio =                   1;
    ODQualityThresholds.ODResults.MinLUPIRatio =                1;
    ODQualityThresholds.EpochAgeTable =                        [0   10
                                                                1   10
                                                                2   5
                                                                3   5
                                                                4   5
                                                                5   5
                                                                6   5
                                                                7   5
                                                                8   3
                                                                9   3
                                                                10  3];
    ODQualityThresholds.EpochAgeDaysToTCA =                     3;
    ODQualityThresholds.Covariance.DefaultSizeFactor =          4.0E+15;
    ODQualityThresholds.Covariance.PortionOfRev =               0.125;
    ODQualityThresholds.CompositeScoreWeightingVector =         [3 2 2 1];
end
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 24-07-2019 |  Initial Development
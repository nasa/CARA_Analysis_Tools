function [DUT1, PMx, PMy] = GetEOPInfo(JDUTC, useDUT1Rate)
% GetEOPInfo - Retrieves Earth orientation parameter information for the
%              Julian Date passed in
%
% Syntax: [DUT1, PMx, PMy] = GetEOPInfo(JDUTC);
%         [DUT1, PMx, PMy] = GetEOPInfo(JDUTC, useDUT1Rate);
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   JDUTC       - Julian Date of the UTC epoch to lookup
%   useDUT1Rate - (Optional) When set to false the DUT1 value is estimated
%                 using linear interpolation from the EOP.mat file. When
%                 set to true, the DUT1 value is estimated using the DUT1
%                 and DUT1Rate values from the time_constants.dat file.
%                 Boolean, defaults to false.
%
% =========================================================================
%
% Output:
%
%   DUT1    - Time difference between UTC and UT1 (universal time) in
%             seconds
%   PMx     - Polar motion x offset in radians
%   PMy     - Polar motion y offset in radians 
%
% =========================================================================
%
% Description:
%
%   This function retrieves the Earth Orientation Parameters (EOP) for the
%   UTC Julian Date passed in. The values are retrieved from the EOP.mat
%   and time_constants.dat files that are available in the top-level
%   DataFiles/TimeCons directory. It is a rewrite/combination of the
%   previous DeltaUT1.m and PolarMotion.m files and adds the functionality
%   of using the time_constants.dat file and DUT1 estimation from DUT1Rate.
%
%   The EOP.mat file is a Matlab formatted file which contains time, DUT1,
%   xp, and yp values throughout time. This format retains all the
%   precision from source data.
%
%   The time_constants.dat file is a fixed-width formatted file which also
%   includes leap seconds and DUT1Rate information. This format loses
%   precision in the DUT1, xp, and yp values; thus it is not the preferred
%   format. However, some historical applications (namely VCM generation)
%   use this format and it is retained to allow replication of those
%   program outputs. Generally, this data will only be used if the date
%   passed in is not covered by the EOP.mat file or if the useUT1Rate
%   option is set to true (which indicates time_constants.dat data is
%   needed).
%
% =========================================================================
%
% Dependencies:
%
%   Data Files: EOP.mat and time_constants.dat
%
% =========================================================================
%
% Initial version: May 2025;  Latest update: May 2025
%
% ----------------- BEGIN CODE -----------------

    %% Check input arguments
    if nargin == 1
        useDUT1Rate = false;
    elseif nargin ~= 2
        error('Incorrect number of arguments passed in!');
    end

    %% Load EOPData from EOP.mat and time_constants.dat
    persistent EOPData;
    
    % Offset between the Julian Date epoch and Matlab's datenum epoch
    jd_datenum_offset = 1721058.5;
    if (isempty(EOPData))
        [p,~,~] = fileparts(mfilename('fullpath'));
        
        % Pull data from both files
        timeConstantsFile = fullfile(p, '../../../DataFiles/TimeCons/time_constants.dat');
        eopFile = fullfile(p, '../../../DataFiles/TimeCons/EOP.mat');

        % Time Constants .dat file contents:
        % yr (2-digit), day of year, date (2-digit year), leap seconds, deltaUT1, deltaUT1Rate, polar x, polar y
        timeConstantsTable = readtable(timeConstantsFile);
        tconsYr = timeConstantsTable.Var1;
        adj1900 = tconsYr > 72;
        tconsYr(adj1900)  = tconsYr(adj1900) + 1900;
        tconsYr(~adj1900) = tconsYr(~adj1900) + 2000;
        tconsDOY = timeConstantsTable.Var2;
        tconsDatenum = datenum(tconsYr,1,tconsDOY); %#ok<DATNM>
        tconsJD = tconsDatenum + jd_datenum_offset;

        % EOPInfo .mat file contents:
        % GMJD (time since 01/05/1941 12:00:00 UTC), polar x, polar y, deltaUT1Rate
        load(eopFile,'EOPInfo');
        eopJD = EOPInfo(:,1) + 2400000.5 + 29999.5;

        % Setup variables to hold combined data
        JD = sort(unique([tconsJD; eopJD]));
        EpochTime = datetime(JD - jd_datenum_offset,'ConvertFrom','datenum');
        tconsMembers = ismember(JD, tconsJD);
        eopMembers = ismember(JD, eopJD);
        numRows = length(JD);
        DeltaUT1          = nan(numRows,1);
        PolarMotionX      = nan(numRows,1);
        PolarMotionY      = nan(numRows,1);
        TconsDeltaUT1     = nan(numRows,1);
        TconsDeltaUT1Rate = nan(numRows,1);

        % Combine the data from the two tables
        eopIndx = 0;
        tconsIndx = 0;
        for i = 1:numRows
            if eopMembers(i)
                eopIndx = eopIndx + 1;
            end
            if tconsMembers(i)
                tconsIndx = tconsIndx + 1;
            end
            if eopMembers(i) && tconsMembers(i)
                % We prefer the EOP.mat information since it retains all
                % digits of precision. However, the tcons data is the
                % format used by ASW, and we should use that data for
                % ASW-specific transformations
                DeltaUT1(i)          = EOPInfo(eopIndx,4);
                PolarMotionX(i)      = EOPInfo(eopIndx,2);
                PolarMotionY(i)      = EOPInfo(eopIndx,3);
                TconsDeltaUT1(i)     = timeConstantsTable.Var5(tconsIndx);
                TconsDeltaUT1Rate(i) = timeConstantsTable.Var6(tconsIndx);
            elseif eopMembers(i)
                DeltaUT1(i)          = EOPInfo(eopIndx,4);
                PolarMotionX(i)      = EOPInfo(eopIndx,2);
                PolarMotionY(i)      = EOPInfo(eopIndx,3);
            elseif tconsMembers(i)
                DeltaUT1(i)          = timeConstantsTable.Var5(tconsIndx);
                PolarMotionX(i)      = timeConstantsTable.Var7(tconsIndx);
                PolarMotionY(i)      = timeConstantsTable.Var8(tconsIndx);
                TconsDeltaUT1(i)     = timeConstantsTable.Var5(tconsIndx);
                TconsDeltaUT1Rate(i) = timeConstantsTable.Var6(tconsIndx);
            end
        end
        EOPData = table(JD,EpochTime,DeltaUT1,PolarMotionX,PolarMotionY,TconsDeltaUT1,TconsDeltaUT1Rate);
    end

    arcsec2rad = (1/3600) * (pi/180);

    % Check for out of range values
    if JDUTC < EOPData.JD(1) || JDUTC > EOPData.JD(end)
        inDate = datetime(JDUTC - jd_datenum_offset,'ConvertFrom','datenum');
        warnTxt = ['Input Julian Date (' num2str(JDUTC) ', ' datestr(inDate,'yyyy-mm-dd HH:MM:SS.FFF') ') is outside of EOP info range:' newline ...
            '  ' datestr(EOPData.EpochTime(1),'yyyy-mm-dd') ' to ' datestr(EOPData.EpochTime(end),'yyyy-mm-dd') newline ...
            'Using closest EOP parameters instead of actual values']; %#ok<*DATST>
        if JDUTC > EOPData.JD(end)
            warnTxt = [warnTxt newline '*** EOP Update is needed! ***'];
        end
        warning(warnTxt);

        if JDUTC < EOPData.JD(1)
            DUT1 = EOPData.DeltaUT1(1);
            PMx  = EOPData.PolarMotionX(1) * arcsec2rad;
            PMy  = EOPData.PolarMotionY(1) * arcsec2rad;
        else
            DUT1 = EOPData.DeltaUT1(end);
            PMx  = EOPData.PolarMotionX(end) * arcsec2rad;
            PMy  = EOPData.PolarMotionY(end) * arcsec2rad;
        end
    else
        % If within the date range, first check for an exact match
        tmpData = EOPData(EOPData.JD == JDUTC,:);
        if ~isempty(tmpData)
            DUT1 = tmpData.DeltaUT1(1);
            PMx  = tmpData.PolarMotionX(1) * arcsec2rad;
            PMy  = tmpData.PolarMotionY(1) * arcsec2rad;
        else
            % No exact match, find indexes bracketing the input date
            Ind = EOPData.JD <= JDUTC;
            idx1 = sum(Ind);
            idx2 = idx1+1;

            % Calculate Polar Motion values
            PMx = LinInterp(EOPData.JD(idx1:idx2),EOPData.PolarMotionX(idx1:idx2),JDUTC) * arcsec2rad;
            PMy = LinInterp(EOPData.JD(idx1:idx2),EOPData.PolarMotionY(idx1:idx2),JDUTC) * arcsec2rad;

            % Estimate UT1 either via interpolation or using the rate
            if ~useDUT1Rate || isnan(EOPData.TconsDeltaUT1Rate(idx1))
                DeltaUT1Vals = EOPData.DeltaUT1(idx1:idx2);
                % Adjust for leap seconds
                if DeltaUT1Vals(2)-DeltaUT1Vals(1) > 0.9
                    % Leap second added
                    DeltaUT1Vals(2) = DeltaUT1Vals(2) - 1;
                elseif DeltaUT1Vals(1)-DeltaUT1Vals(2) > 0.9
                    % Leap second subtracted
                    DeltaUT1Vals(2) = DeltaUT1Vals(2) + 1;
                end
                DUT1 = LinInterp(EOPData.JD(idx1:idx2),DeltaUT1Vals,JDUTC);
            else
                % This method is specific to ASW, so use the Tcons data
                % instead
                % Time difference in days
                timeDiff = JDUTC - EOPData.JD(idx1);
                % Converts from milliseconds/day to seconds/day
                ratePerDay = EOPData.TconsDeltaUT1Rate(idx1) / 1000;
                DUT1 = EOPData.TconsDeltaUT1(idx1) + (timeDiff * ratePerDay);
            end
        end
    end
end

function [Y0] = LinInterp(X,Y,X0)

    % Compute slope
    m  = (Y(2)-Y(1)) / (X(2)-X(1));
    
    % Linear Interpolation
    Y0 = Y(1) + m * (X0 - X(1));

end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% L. Baars       | 05-13-2025 | Combined the DeltaUT1.m and PolarMotion.m
%                               files into this function since most of the
%                               calculations are equivalent. Added the
%                               functionality to also read in
%                               time_constants.dat and provide a DUT1
%                               estimation option using the DUT1 and
%                               DUT1Rate values from time_constants.dat.

% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

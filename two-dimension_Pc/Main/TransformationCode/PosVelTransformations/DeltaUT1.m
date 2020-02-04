function [DUT1] = DeltaUT1(JDUTC)
%
% DeltaUT1 - This function computes the difference between UTC and UT1 at a
%            specific epoch. Necessary data for the computatation provided
%            in the EOP.mat file. If the epoch doesn't fall on the exact date
%            in the EOP file then the appropriate value is determined by
%            linearly interpolating two adjacent values.
%
% Syntax:   [DUT1] = DeltaUT1(JDUTC)
%
% Inputs:
%   JDUTC  - Julian Date of the epoch (UTC) [Scalar]
%
% Outputs:
%   DUT1   - Delta time between UTC and UT1 [sec]
%            Never exceeds 0.9 seconds
%
% Examples/Validation Cases: 
%
% Other m-files required: None
% Subfunctions: LinInterp
% MAT-files required: EOP.mat
% Global variables: EOPInfo
%
% See also: None
%
% July 2015; Last revision: 15-May-2016
%
% ----------------- BEGIN CODE -----------------
    
    % EOP data
    global EOPInfo;
     
    if (isempty(EOPInfo))
        % Get EOP info from a .mat file
        load([cd '\EOP.mat']);
    end

    % Modified Julian Date (GSFC)
    MJD_GSFC    = EOPInfo(:,1);
    MJDUTC_GSFC = JDUTC - 2400000.5 - 29999.5;
    
    % Check if JDUTC is out of range of applicable values
    if (MJDUTC_GSFC < MJD_GSFC(1) || MJDUTC_GSFC > MJD_GSFC(end))
        
        fprintf('EpochUTC is out of range. Consider updating EOP.mat file.\n');
        fprintf('Closest EOP value used for the calculation.\n');
        
        if (MJDUTC_GSFC < MJD_GSFC(1))
        
            DUT1 = EOPInfo(1,4);
            
            return
 
        elseif (MJDUTC_GSFC > MJD_GSFC(end))
            
            DUT1 = EOPInfo(end,4);

            return
            
        end
                       
     end
    
    % Find index if epoch matches any of the EOP dates perfectly
    Id = find(MJD_GSFC == MJDUTC_GSFC,1);
    
    % If yes
    if (~isempty(Id))
        
        % No interpolation necessary here
        DUT1 = EOPInfo(Id,4);
        
    % If no, then interpolate between two values (linearly)
    else
    
        Id = MJD_GSFC <= MJDUTC_GSFC;  
    
        x1 = EOPInfo(sum(Id),1);
        x2 = EOPInfo(sum(Id)+1,1);
        
        y1 = EOPInfo(sum(Id),4);
        y2 = EOPInfo(sum(Id)+1,4);
        
        % Subtract 1 sec UT1 if y2 is a leap sec day and y1 is before the
        % leapsecond
        if y2 > 0 && y1 < 0
            y2 = y2 - 1;
        end
     
        % Linear interpolation
        [DUT1] = LinInterp([x1,x2],[y1,y2],MJDUTC_GSFC);
        
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
% D. Plakalovic  | 07-13-2015 |  Re-coded this function from the original 
%                                version (Feb 2013). Inserted additional
%                                functionality.
% L. Johnson     | 05-15-2016 |  Modified code to account for records over
%                                a leapsecond, i.e. the days before the day
%                                of a leadpsecond
%             
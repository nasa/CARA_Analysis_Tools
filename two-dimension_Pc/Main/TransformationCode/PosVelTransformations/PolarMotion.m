function [xp,yp] = PolarMotion(JDUTC)
%
% PolarMotion - This function computes the polar motion parameters at a
%               specific epoch. Necessary data for the computatation provided
%               in the EOP.mat file. If the epoch doesn't fall on the exact date
%           	in the EOP file then the appropriate value is determined by
%               linearly interpolating two adjacent values.
%
% Syntax:      [xp,yp] = PolarMotion(JDUTC)
%
% Inputs:
%   JDUTC     - Julian Date of the epoch (UTC) [Scalar]
%
% Outputs:
%   xp        - Angle describing the Earth's spin axis location as it is 
%               not fixed in the Earth [radians]
%   yp        - Angle describing the Earth's spin axis location as it is 
%               not fixed in the Earth [radians]
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
% July 2015; Last revision: 13-Jul-2015
%
% ----------------- BEGIN CODE -----------------

    % EOP data
    global EOPInfo;
     
    if (isempty(EOPInfo))
        % Get EOP info from a .mat file
        load('EOP.mat');
    end
    
    % Modified Julian Date (GSFC)
    MJD_GSFC    = EOPInfo(:,1);
    MJDUTC_GSFC = JDUTC - 2400000.5 - 29999.5;
    
    % Check if JDUTC is out of range of applicable values
    if (MJDUTC_GSFC < MJD_GSFC(1) || MJDUTC_GSFC > MJD_GSFC(end))
        
        fprintf('EpochUTC is out of range. Consider updating EOP.mat file.\n');
        fprintf('Closest EOP value used for the calculation.\n');
        
        if (MJDUTC_GSFC < MJD_GSFC(1))
        
            xp = EOPInfo(1,2)*(1/3600)*(pi/180);
            yp = EOPInfo(1,3)*(1/3600)*(pi/180);
            
            return
 
        elseif (MJDUTC_GSFC > MJD_GSFC(end))
            
            xp = EOPInfo(end,2)*(1/3600)*(pi/180);
            yp = EOPInfo(end,3)*(1/3600)*(pi/180);
            
            return
            
        end
                       
     end
        
    % Find index if epoch matches any of the EOP dates perfectly
    Id = find(MJD_GSFC == MJDUTC_GSFC,1);
    
    % If yes
    if (~isempty(Id))
        
        % No interpolation necessary here
        xp = EOPInfo(Id,2)*(1/3600)*(pi/180);
        yp = EOPInfo(Id,3)*(1/3600)*(pi/180);
        
    % If no, then interpolate between two values (linearly)
    else
    
        Id   = MJD_GSFC <= MJDUTC_GSFC;
    
        x1   = EOPInfo(sum(Id),1);
        x2   = EOPInfo(sum(Id)+1,1);
        
        xp1  = EOPInfo(sum(Id),2)  *(1/3600)*(pi/180);
        xp2  = EOPInfo(sum(Id)+1,2)*(1/3600)*(pi/180);
        
        yp1  = EOPInfo(sum(Id),3)  *(1/3600)*(pi/180);
        yp2  = EOPInfo(sum(Id)+1,3)*(1/3600)*(pi/180);
        
        [xp] = LinInterp([x1,x2],[xp1,xp2],MJDUTC_GSFC);
        [yp] = LinInterp([x1,x2],[yp1,yp2],MJDUTC_GSFC);
        
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
% L. Johnson     | 02-23-2016 |  Modified code to always return angles in
%                                radians
%           
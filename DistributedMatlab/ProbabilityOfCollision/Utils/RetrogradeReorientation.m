function [r1,v1,C1,r2,v2,C2,out] = RetrogradeReorientation(r1,v1,C1,r2,v2,C2,params)

% Reorient primary and secondary states to prevent either from
% being processed as a retrograde Keplerian orbit (i.e., an orbit 
% with 180deg inclination). This entails rotating both states to use
% a frame with a different set of inertial axes. Also reorient the
% covariance matrices.
%
% Input:
%
%  (r1,v1) = Primary   position velocity  [1x3 or 3x1 vectors]
%
%  (r2,v2) = Secondary position velocity  [1x3 or 3x1 vectors]
%
%  (C1,C2) = Primary and secondary covariances [6x6 matrices]
%
%  params.RetrogradeReorientation = processing mode
%   0 => No retrograde orbit reorientation
%   1 => If either orbit is retrograde, reorient the ref. frame axes
%   2 => Always reorient the ref. frame axes (testing mode)
%   3 => Reorient axes to force primary to be retrograde (testing mode)
%
% Output:
%
%  (r1,v1) = Primary   position velocity  [1x3 or 3x1 vectors]
%
%  (r2,v2) = Secondary position velocity  [1x3 or 3x1 vectors]
%
%  (C1,C2) = Primary and secondary covariances [6x6 matrices]
%
%  out = Output strcuture:
%  out.Reoriented = Binary flag to indicate if frame reorientation
%                   was performed
%  out.Retro1 = Primary   orbit retrograde flag, if calculated
%  out.Retro2 = Secondary orbit retrograde flag, if calculated
%

% Set defaults
if nargin < 7; params = []; end
params = set_default_param(params,'RetrogradeReorientation',1);

% Process different modes
if params.RetrogradeReorientation == 0
    
    % No retrogade or reorientation processing required
    out.Reoriented = false;
    
else
    
    % Check the sizes of input vectors
    sz = size(r1);
    if ~isequal(sz,size(r2)) || ~isequal(sz,size(v1)) || ~isequal(sz,size(v2))
        error('Input vectors must be same dimension');
    else
        % Set InputRowVectors flag, and transform to column vectors if
        % required
        if isequal(sz,[1 3])
            InputRowVectors = true;
            r1 = reshape(r1,[3 1]);
            v1 = reshape(v1,[3 1]);
            r2 = reshape(r2,[3 1]);
            v2 = reshape(v2,[3 1]);
        elseif isequal(sz,[3 1])
            InputRowVectors = false;
        else
            error('Invalid dimension for input vectors');
        end
    end

    % Process different retrograde reorientation modes
    if params.RetrogradeReorientation == 1 || ...
       params.RetrogradeReorientation == 2
        
        % Check if frame reorientation is required due to orbit being
        % within a small tolerance of retrograde
        
        % Set default tolerance
        params = set_default_param(params,'Eps',1e-6);
        
        % Check if primary orbit is near retrograde
        out.Retro1 = CheckForRetrograde(r1,v1,params.Eps);
        
        % Check if secondary orbit is near retrograde
        out.Retro2 = CheckForRetrograde(r2,v2,params.Eps);
        
        % Adjust reference frame, attempting to eliminate retrograde orbits
        if out.Retro1 || out.Retro2 || params.RetrogradeReorientation == 2
            
            % Set default verbose flag
            params = set_default_param(params,'verbose',true);
            if params.verbose
                disp(['*** Reorienting frame. Retrograde status:' ...
                    ' Primary = ' num2str(out.Retro1) ...
                    ' Secondary = ' num2str(out.Retro2)]);
            end
        
            % Axes for new reference frame (corresponds to Euler angles of
            % phi = 90deg, theta = 90deg, psi = 0).
            Xhat =  [0; 1; 0];
            Yhat =  [0; 0; 1];
            Zhat =  [1; 0; 0];
            
            % Rotation matrix for reorioented ref. frame
            M3 = [Xhat Yhat Zhat]';
            
            % Rotate states to use reoriented frame
            r1New = M3 * r1; v1New = M3 * v1;
            r2New = M3 * r2; v2New = M3 * v2;
            
            % Check that reorientation did not result in retrograde orbits
            Retro1New = CheckForRetrograde(r1New,v1New,params.Eps);
            Retro2New = CheckForRetrograde(r2New,v2New,params.Eps);
            
            % If retrograde orbits persist (should happen very, very
            % rarely), then iterate to find a reorientation that eliminates
            % the retrograde status
            if Retro1New || Retro2New
                
                % Set default max. number of reorientation iterations
                params = set_default_param(params,'MaxIter',5);
                
                % Iterate to find a reorientation 
                iterating = true; iter = 0; converged = false;
                halfpi = pi/2; phi = halfpi; tht = halfpi;
                while iterating
                    disp([' ** Reorientation iteration ' num2str(iter) ...
                        ' for phi = ' num2str(phi*180/pi) ...
                        ' & theta = ' num2str(tht*180/pi) ...
                        ' Retrograde status:' ...
                        ' Primary = ' num2str(Retro1New) ...
                        ' Secondary = ' num2str(Retro2New)]);
                    % Randomly generate new phi and tht Euler angles
                    phi = halfpi*(rand+1.0); cp=cos(phi); sp=sin(phi);
                    tht = halfpi*(rand+0.5); ct=cos(tht); st=sin(tht);
                    % Calculate new reorientation matrix
                    C = [1 0 0; 0 ct st; 0 -st ct];
                    D = [cp sp 0; -sp cp 0; 0 0 1];
                    M3 = C*D;
                    % Check that reorientation did not produce retrogrades
                    Retro1New = CheckForRetrograde(r1New,v1New,params.Eps);
                    Retro2New = CheckForRetrograde(r2New,v2New,params.Eps);
                    if ~Retro1New && ~Retro2New
                        converged = true;
                        iterating = false; 
                    else
                        iter = iter+1;
                        iterating = iter < MaxIter;
                    end
                end
                
            else
                
                % First reorientation converged to non-retrograde status
                converged = true;
                
            end
            
            % If converged, adopt new orientation
            if converged
                % Output new vectors
                r1 = r1New; v1 = v1New; r2 = r2New; v2 = v2New;
                % Output rotation matrices
                out.M3 = M3; % 3x3 rotation matrix
                Z3 = zeros(3,3); out.M6 = [M3 Z3; Z3 M3]; % 6x6 rot. matrix
                % Output new covariances
                C1 = out.M6 * C1 * out.M6'; C2 = out.M6 * C2 * out.M6';
                % Set the frame reorientation flag
                out.Reoriented = true;
            else
                error('No retrograde reorientation found');
            end
            
        else
            
            % No reorientation required
            out.Reoriented = false;
            
        end
        
    elseif params.RetrogradeReorientation == 3
        
        % Special testing mode to reorient such that primary orbit is
        % perfectly retrograde
        out.Reoriented = true;

        % Set default verbose flag
        params = set_default_param(params,'verbose',true);
        if params.verbose
            warning('Reorienting frame to make primary orbit retrograde');
        end

        % Force primary to have a perfectly retrograde orbit (for testing)
        hvec1 = cross(r1,v1); hmag1 = norm(hvec1);
        Zhat = -hvec1/hmag1;
        Xhat = r1-Zhat*(r1'*Zhat); Xhat = Xhat/norm(Xhat);
        Yhat = cross(Zhat,Xhat);
        M3 = [Xhat Yhat Zhat]'; Z3 = zeros(3,3); M6 = [M3 Z3; Z3 M3];

        % Rotate states and covariances
        r1 = M3 * r1; v1 = M3 * v1; C1 = M6 * C1 * M6';
        r2 = M3 * r2; v2 = M3 * v2; C2 = M6 * C2 * M6';
        
    else

        error('Invalid RetrogradeReorientation parameter');

    end
    
    % Convert output to match input row vectors, if required
    if InputRowVectors
        r1 = reshape(r1,[1 3]);
        v1 = reshape(v1,[1 3]);
        r2 = reshape(r2,[1 3]);
        v2 = reshape(v2,[1 3]);
    end    
    
end

return
end

% =========================================================================

function Retro = CheckForRetrograde(r,v,Eps)
% Check if an orbit is near retrograde
h = cross(r,v);
Retro = 1+h(3)/sqrt(h'*h) < Eps;
return
end


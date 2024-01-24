function [cov] = CheckAndResizeCov(numR, cov)
% CheckAndResizeCov - Resizes the covariance passed in into an n x 9
%                     matrix to be used in vectorized covariance
%                     processing.
%
% Description:
%
%   Reformats the input covariance into an n x 9 matrix representing the
%   3x3 position covariance in the following format:
%     cov = [C(1,1) C(1,2) C(1,3) C(2,1) C(2,2) C(2,3) C(3,1) C(3,2) C(3,3)]
%   Where "n" is the number of rows passed in.
%
%   Depending on the input format, this function will respond in several
%   different ways to reformat the cov matrix:
%   - If the input is an nx9, then this function verifies that n matches
%     numR.
%   - If the input is a 3x3, then this function will reformat the matrix
%     into a 1x9 and then repeat this vector numR times to create n rows.
%   - If the input is a 6x6, then this function will take the upper left
%     3x3 component and reformat it into a 1x9. Then it will repeat this
%     vector numR times to create n rows.
%   - If the input is a 3x3xn, then this function will verify that n
%     matches numR and will reformat each 3x3 into a corresponding 1x9.
%   - If the input is a 6x6xn, then this function will verify that n
%     matches and will take the upper left 3x3 component of each 6x6 into a
%     corresponding 1x9.
%   - Any other input is considered an error.

    covSize = size(cov);
    % If a 2D covariance matrix was passed in
    if size(covSize,2) == 2
        if covSize(2) ~= 9 && (covSize(1) ~= 3 || covSize(2) ~= 3) && ...
                              (covSize(1) ~= 6 || covSize(2) ~= 6)
            error('CheckAndResizeCov:Invalid2DCov','2D Covariance matrix must have 9 columns or be a 3x3 or 6x6 matrix!');
        end
        % Resize down to a 3x3 if a 6x6 was passed in
        if covSize(1) == 6 && covSize(2) == 6
            cov = cov(1:3,1:3);
            covSize = size(cov);
        end
        if covSize(1) == 1
            cov = repmat(cov,numR,1);
        elseif covSize(1) == 3 && covSize(2) == 3
            cov = reshape(permute(cov,[2 1]),1,9);
            cov = repmat(cov,numR,1);
        elseif covSize(1) ~= numR
            error('CheckAndResizeCov:RowCountMismatch2D','2D Covariance cannot be resized to match r matrix');
        end
    % If a 3D covariance matrix was passed in
    elseif size(covSize,2) == 3
        % The 3D matrix should have the dimension 3x3xnumR or 6x6xnumR
        if (covSize(1) ~= 3 || covSize(2) ~= 3 || covSize(3) ~= numR) && ...
           (covSize(1) ~= 6 || covSize(2) ~= 6 || covSize(3) ~= numR)
            error('CheckAndResizeCov:Invalid3DCov','3D covariance matrix must be of size 3x3xnumR or 6x6xnumR');
        end
        % Resize down to a 3x3xnumR if a 6x6xnumR was passed in
        if covSize(1) == 6 && covSize(2) == 6
            cov = cov(1:3,1:3,:);
        end
        cov = reshape(permute(cov,[3 2 1]),numR,9);
    else
        error('CheckAndResizeCov:InvalidCov','Improperly sized covariance was detected');
    end
end
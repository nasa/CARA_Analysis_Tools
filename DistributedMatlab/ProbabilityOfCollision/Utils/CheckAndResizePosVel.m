% checks the sizes of the matrices and resizes v if needed
function [numR, v] = CheckAndResizePosVel(r, v)
    [numR, rColumns] = size(r);
    if rColumns ~= 3
        error('r matrix must have 3 columns!');
    end
    [numV, vColumns] = size(v);
    if vColumns ~= 3
        error('v matrix must have 3 columns!');
    end
    if numV ~= numR
        if numV == 1
            v = repmat(v,numR,1);
        else
            error('v matrix cannot be resized to match r matrix');
        end
    end
end
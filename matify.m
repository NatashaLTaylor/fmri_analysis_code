function out_mat = matify(vec,nROI)

%MATIFY       Turn a vectorized version of the upper triangle of matrix
%             (minus the diagonal) into a symmetrical matrix
%
%   Inputs:
%       vec,
%           an Nx1 vector, where N is the number of unique values from the
%           upper triangle (minus the diagonal) of an adjacency matrix
%       nROI, number of regions for matrix
%
%   Outputs:
%       out_mat
%           a symmetric adjacency matrix


    % create a lookup table of all possible matrix sizes (up to N = 400)
    lookup_table = transpose(1:1:nROI);
    for nn = 1:nROI
        lookup_table(nn,2) = (nn*nn-nn)/2;
    end
    
    % find out how many nodes should be in the matrix
    size1 = size(vec,1);
    vec_idx = lookup_table(:,2)==size1;
    nNodes = lookup_table(vec_idx,1);
    out_mat = zeros(nNodes);
    
    % fill in the matrix
    tally = 0;
    for xx = 1:nNodes-1
        temp = size(out_mat(xx,xx+1:end),2);
        out_mat(xx,xx+1:end) = vec(tally+1:temp+tally);
        tally = temp+tally;
    end

    % make the bottom left look like the top right
    out_mat = triu(out_mat) + triu(out_mat)';

end
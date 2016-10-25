function [D, ni] = find_nn(X, k, method)
%FIND_NN Finds k nearest neigbors for all datapoints in the dataset
%
%	[D, ni] = find_nn(X, k, method)
%
% Finds the nearest neighbors for all datapoints in the dataset X.
% In X, rows correspond to the observations and columns to the
% dimensions. The value of k is the number of neighbors that is
% stored. The function returns a sparse distance matrix D, in which
% only the distances to the nearest neighbors are stored. For
% equal datapoints, the distance is set to a tolerance value.
% The method is relatively slow, but has a memory requirement of O(nk).
%
% Method:
%       neighbors       find the k nearest neighbors
%       region          find the nearest neighbors in a region at a
%                       distance k

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('method', 'var') || ~ismember(method,{'neighbors', 'region'})
        method = 'neighbors';
    end
    
    if ~exist('k', 'var') || isempty(k)
        k = 12;
    end
    
    switch method
        case 'neighbors'
            % Perform adaptive neighborhood selection if desired
            if ischar(k)
                [D, max_k] = find_nn_adaptive(X);
                if (nargout > 1)
                    ni = zeros(size(X, 1), max_k);
                    for i=1:size(X, 1)
                        tmp = find(D(i,:) ~= 0);
                        ni(i,1:length(tmp)) = tmp;
                    end
                end

            % Perform normal neighborhood selection
            else
                % Compute distances in batches
                n = size(X, 1);
                sum_X = sum(X .^ 2, 2);
                batch_size = round(2e7 ./ n);
                D = zeros(n, k);
                ni = zeros(n, k);
                for i=1:batch_size:n
                    batch_ind = i:min(i + batch_size - 1, n);
                    DD = bsxfun(@plus, sum_X', bsxfun(@plus, sum_X(batch_ind), ...
                                                           -2 * (X(batch_ind,:) * X')));
                    [DD, ind] = sort(abs(DD), 2, 'ascend');
                    D(batch_ind,:) = sqrt(DD(:,2:k + 1));
                    ni(batch_ind,:) = ind(:,2:k + 1);
                end
                D(D == 0) = 1e-9;
                % D = sparse(repmat(1:n, [1 k])', ni(:), D(:), n, n);
                % D = sparse([repmat(1:n, [1 k])'; ni(:)], [ni(:); repmat(1:n, [1 k])'], [D(:); D(:)], n, n);
                % fix the sparse bug: doubling values for mutual nearest neighbor
                % nodes
                Dout = zeros(n,n);
                idx = repmat(1:n, [1 k])';
                Dout(sub2ind([n,n],idx,ni(:))) = D;
                Dout(sub2ind([n,n],ni(:),idx)) = D;
                D = sparse(Dout);
            end
        case 'region'
            % Compute all distances: this is memory demanding but faster
            % than squareform(pdist(X))
            DD = sum(X .^ 2);
            D = real(sqrt(bsxfun(@plus, DD.', DD) - (2 * (X.' * X))));
            D(1:size(X, 1)+1:end) = 0;
            D(D > K) = 0;
            if (nargout > 1)
                max_k = max(sum(D>0, 2));
                ni = zeros(size(X, 1), max_k);
                for i=1:size(X, 1)
                    tmp = find(D(i,:) ~= 0);
                    ni(i,1:length(tmp)) = tmp;
                end
            end
            D = sparse(D);
    end
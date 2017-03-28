function [mappedX, mapping] = isomap(X, no_dims, k, nn_method)
%ISOMAP Runs the Isomap algorithm
%
%   [mappedX, mapping] = isomap(X, no_dims, k, nn_method); 
%
% The functions runs the Isomap algorithm on dataset X to reduce the
% dimensionality of the dataset to no_dims. The number of neighbors used in
% the compuations is set by k (default = 12). This implementation does not
% use the Landmark-Isomap algorithm.
% nn_method is the nearest neighbor method: region for epsilon NN and
% neighbors for traditional NN (default).
%
% If the neighborhood graph that is constructed is not completely
% connected, only the largest connected component is embedded. The indices
% of this component are returned in mapping.conn_comp.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('no_dims', 'var')
        no_dims = 2;
    end
    if ~exist('k', 'var')
        k = 12;
    end
    if ~exist('nn_method', 'var')  || ~ismember(nn_method,{'neighbors', 'region'})
        nn_method = 'neighbors';
    end
    
    % Construct neighborhood graph
    disp('Constructing neighborhood graph...');    
  D = find_nn(X, k, 'region');
 
    
    % Select largest connected component
    blocks = components(D)';
    count = zeros(1, max(blocks));
    for i=1:max(blocks)
        count(i) = length(find(blocks == i));
    end
    [~, block_no] = max(count);
    conn_comp = (blocks == block_no);    
    D = D(conn_comp, conn_comp);
    mapping.D = D;
    n = size(D, 1);
    % Compute shortest paths
    disp('Computing shortest paths...');
    D = dijkstra(D, 1:n);
    mapping.DD = D;
    % Performing MDS using eigenvector implementation
    disp('Constructing low-dimensional embedding...');
    D = D .^ 2;
    M = -.5 .* (bsxfun(@minus, bsxfun(@minus, D, sum(D, 1)' ./ n), sum(D, 1) ./ n) + sum(D(:)) ./ (n .^ 2));
    M(isnan(M)) = 0;
    M(isinf(M)) = 0;
    M = (M+M')/2;
    dev = gpuDevice;
    useGpu = ~isempty(dev) && (dev.AvailableMemory > 5e9); % 5GB
    if (useGpu)
        [vec, val] = eig(gpuArray(single(M)));
        vec = gather(vec);
        val = gather(val);
    else
        [vec, val] = eig(M);
    end
    if (size(vec, 2) < no_dims)
		no_dims = size(vec, 2);
		warning(['Target dimensionality reduced to ' num2str(no_dims) '...']);
    end
    % Computing final embedding
    [val, ind] = sort(real(diag(val)), 'descend'); 
    vec = vec(:,ind(1:no_dims));
    val = val(1:no_dims);
    mappedX = real(bsxfun(@times, vec, sqrt(val)'));
    
    % Store data for out-of-sample extension
    mapping.conn_comp = conn_comp;
    mapping.k = k;
    mapping.X = X(conn_comp,:);
    mapping.vec = vec;
    mapping.val = val;
    mapping.no_dims = no_dims;
    
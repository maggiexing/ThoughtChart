
function R = isomapResidualVariance(mapping)
% Compute the residual variance for the embedding based on the output of the 
% isomap method found in the Matlab Toolbox for Dimensionality Reduction.

% R = isomapResidualVariance(mapping)
% Input:
% mapping is a struct made of
%             X: original matrix (N X M) N observations and M features
%             D: nearest neighbor matrix (N X N)
%            DD: Dijkstra's shortest path based on D (N X N)
%     conn_comp: connected components (N X 1)
%             k: nearest neighbors
%           vec: eigen vector (N X no_dims)
%           val: eigen value (no_dims X 1)
%       no_dims: number of dimensions X should be reduced to
%          name: 'Isomap'
% Output: 
%    R = residual variances for embeddings in Y


%%
N = size(mapping.X,1);
D = reshape(mapping.D,N^2,1);
dims = mapping.no_dims;
R = zeros(dims, 1);
for di = 1:dims
     if (di<=N)
         Y = (real(mapping.vec(:,1:di)).*(ones(N,1)*sqrt(mapping.val(1:di))'))';
         r2 = 1-corrcoef(reshape(real(L2_distance(Y, Y, 1)),N^2,1),D).^2; 
         R(di) = r2(2,1); 
     end
end

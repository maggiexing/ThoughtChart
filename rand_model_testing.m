% %% This is the code for checking the spiral shapeCompareH1 for null data set
%  %gpuDevice;
% %   % Load the data into one matrix
gpuDevice;
%% creat 10000 matrix goes from all zeros to all 1

val=0:1/9999:1;
DyMat=zeros(34,34,length(val));
for i=1:length(val)
    mat(1:34,1:34) = val(i);
    DyMat(:,:,i)=mat;
end
%add noise
DyMat=DyMat+normrnd(0.00001,0.0025,34,34,length(val));
for j=1:34
    for k=1:34
        DyMatAll(j,k,:)= DyMat(k,j,:);
        DyMatAll(j,j,:)=0;
    end
end
%% Load the data into one matrix
Sample = length(val);
CnctDim = 34;
%% Run dissimilarity here
% compute Euclidean distance from each pair of connectomes
X = reshape(DyMatAll,[CnctDim*CnctDim,length(val)]);
D = sum(X .^ 2);
Space = real(sqrt(bsxfun(@plus, D.', D) - (2 * (X.' * X))));
Space(1:length(Space)+1:end) = 0;
clear X D % clear memory from large arrays
% Run NDR here
Edim=10;

[test, dumpAll]=compute_mapping(Space,'Isomap', 10,40);
figure;
subplot(1,3,3)
plot3(test(:,1), test(:,2), test(:,3), '.');
grid on
subplot(1,3,1)
imagesc(Space)
axis equal
axis image
subplot(1,3,2)
hist( DyMatAll(:),500)
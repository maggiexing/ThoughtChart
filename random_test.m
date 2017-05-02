gpuDevice;
clear all
mat= normrnd(0.5,0.25,34,34,10000);
for i=1:34
    for j=1:34
        mat(i,j,:)=mat(j,i,:);
    end
end
mat(1:35:end) = 0;
%% Run dissimilarity here
% compute Euclidean distance from each pair of connectomes
X = reshape(mat,[34*34,10000]);
% L2
D = sum(X .^ 2);
DyDistAll = real(sqrt(bsxfun(@plus, D.', D) - (2 * (X.' * X))));
DyDistAll(1:length(DyDistAll)+1:end) = 0;
clear X D % clear memory from large arrays
%% Run NDR here
Edim=10;
[IsomapXYZ, dumpAll]=compute_mapping(DyDistAll ,'Isomap', Edim,60,'neighbors');
%% plots
IsomapAll=zeros(length(dumpAll.conn_comp),3);
IsomapAll(dumpAll.conn_comp,:)=IsomapXYZ(:,1:3);
figure;
subplot(1,3,3)
plot3(IsomapAll(:,1), IsomapAll(:,2), IsomapAll(:,3), '.');
grid on
subplot(1,3,1)
imagesc(DyDistAll)
axis equal
axis image
subplot(1,3,2)
hist(mat(:),500)
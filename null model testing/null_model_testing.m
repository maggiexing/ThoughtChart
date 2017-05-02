% %% This is the code for checking the spiral shapeCompareH1 for null data set
%  %gpuDevice;
% %   % Load the data into one matrix
gpuDevice;
load WPLI_Neutral_DZ30.mat
load WPLI_Neutral_HC1.mat
%% pick the head matrix and tail matrix from the group
% the fourth neutral in first HC subject
matrix1=mean(abs(Neutral_Control(:,:,4:7,4)),3);
% the 12th DZ subject's neutral at 11th points
matrix2=mean(abs(Neutral_Disease(:,:,4:7,11)),3);
sz = [34 34];       
%% Interpolation 
Q = cell(10000,1);
Q{1}=matrix1;
Q{10000}=matrix2;
% indices of missing matrices
idx = ~cellfun(@isempty,Q);
x = 1:numel(Q);

% merge cells into a multidimensional matrix, call INTERP1, 
QQ = Q(idx);
QQ = permute(cat(3,QQ{:}), [3 1 2]);
QQ = interp1(x(idx), QQ, x);            % one call to interpolation function
QQ = reshape(num2cell(permute(QQ, [2 3 1]), [1 2]), 1,[]);

%% Load the data into one matrix
SampleSizeHC = 10000;
NEEGPoints = 1;
CnctDim = 34;
NPointsHC = NEEGPoints*SampleSizeHC;

TotalNPoints = NPointsHC;
%control group
DyMatAll = zeros(CnctDim, CnctDim, TotalNPoints);
for subjId = 1: SampleSizeHC
    DyMatAll(:,:, (1+NEEGPoints*(subjId-1)):(NEEGPoints*(subjId))) = abs(QQ{1,subjId});
end

%% add noise 
noise_column=randi([1 34],4,1);
noise_row=randi([1 34 ],1,4);
DyMatAll(noise_column,:,:)=rand(4,34,10000);
DyMatAll(:,noise_row,:)=rand(34,4,10000);
DyMatAll(1:NEEGPoints+1:end) = 0;
%% Run dissimilarity here
% compute Euclidean distance from each pair of connectomes
X = reshape(DyMatAll,[CnctDim*CnctDim,TotalNPoints]);
D = sum(X .^ 2);
Space = real(sqrt(bsxfun(@plus, D.', D) - (2 * (X.' * X))));
Space(1:length(Space)+1:end) = 0;
clear X D % clear memory from large arrays
% Run NDR here
Edim=10;
[test, dumpAll]=compute_mapping(Space,'Isomap', 10,40);
plot3(test(:,3), test(:,4), test(:,5), 'o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',2);
grid on;
% % This is the main processing code for dissimilarity based Thought Chart
% for subject conducting ERT task, there are two groups HC vs SAD subjects
% in three conditions: neutral, maintain and reappraise. Input are dynmaic
% connectivity matrix each in size 34*34*50*130, and foi averaged matrix are
% taken for dissimilarity matrix construction.
% Noted that Neutral is alway marked in green
% maintain is in magenta and reappraise is in blue
% script write in 05/29/2016 by Mengqi Xing
%
% init the gpu the first time, this may take some time if running for the
% first time: Driver is compiling binaries for matlab if never done before
% define the following in enviornment variables:
% CUDA_CACHE_MAXSIZE 2147483648
% CUDA_CACHE_DISABLE 0
gpuDevice;
clear all

% loading and finding the foi average for 2 groups
ConnectC={'1','6','9','13','17','19','22','35','39','44','47','49','51','52','53','54','57'};
ConnectP={'8','17','21','23','30','37','40','50','118','124','127','138','140','143','149','153','167'};
foi=4:7;%define the frequency of interest

dataPath = 'F:\Data\Resting&ERT\';
% control
for i=1:length(ConnectC)
    RestC=[dataPath 'WPLI_Reappraise_HC' ConnectC{i} '.mat'];
    R=cell2mat(struct2cell(load(RestC)));
    Alpha_R(:,:,:)=abs(mean(R(:,:,foi,1:130),3));
    ControlA_R{i}=Alpha_R;
    
    RestC=[dataPath 'WPLI_Neutral_HC' ConnectC{i} '.mat'];
    N=cell2mat(struct2cell(load(RestC)));
    Alpha_N(:,:,:)=abs(mean(N(:,:,foi,1:130),3));
    ControlA_N{i}=Alpha_N;
    
    RestC=[dataPath 'WPLI_Maintain_HC' ConnectC{i} '.mat'];
    M=cell2mat(struct2cell(load(RestC)));
    Alpha_M(:,:,:)=abs(mean(M(:,:,foi,1:130),3));
    ControlA_M{i}=Alpha_M;
    
    RestC=[dataPath 'CSAD' ConnectC{i} 'F_CW.mat'];
    R=cell2mat(struct2cell(load(RestC)));
    Alpha_R1(:,:,:)=abs(mean(R(:,:,foi,1:70),3));
    ControlA_CW{i}=Alpha_R1;
    
    RestC=[dataPath 'CSAD' ConnectC{i} 'A_R.mat'];
    R=cell2mat(struct2cell(load(RestC)));
    Alpha_R2(:,:,:)=abs(mean(R(:,:,foi,1:110),3));
    ControlA_AR{i}=Alpha_R2;
end

% disease
for i=1:length(ConnectP)
    RestP=[dataPath 'WPLI_Reappraise_DZ' ConnectP{i} '.mat'];
    R=cell2mat(struct2cell(load(RestP)));
    Alpha_R(:,:,:)=mean(R(:,:,foi,1:130),3);
    DiseaseA_R{i}=Alpha_R;
    
    RestP=[dataPath 'WPLI_Neutral_DZ' ConnectP{i} '.mat'];
    N=cell2mat(struct2cell(load(RestP)));
    Alpha_N(:,:,:)=mean(N(:,:,foi,1:130),3);
    DiseaseA_N{i}=Alpha_N;
    
    RestP=[dataPath 'WPLI_Maintain_DZ' ConnectP{i} '.mat'];
    M=cell2mat(struct2cell(load(RestP)));
    Alpha_M(:,:,:)=mean(M(:,:,foi,1:130),3);
    DiseaseA_M{i}=Alpha_M;
    
    RestC=[dataPath 'PSAD' ConnectP{i} 'F_CW.mat'];
    R=cell2mat(struct2cell(load(RestC)));
    Alpha_R1(:,:,:)=abs(mean(R(:,:,foi,1:70),3));
    DiseaseA_CW{i}=Alpha_R1;
    
    RestC=[dataPath 'PSAD' ConnectP{i} 'A_R.mat'];
    R=cell2mat(struct2cell(load(RestC)));
    Alpha_R2(:,:,:)=abs(mean(R(:,:,foi,1:110),3));
    DiseaseA_AR{i}=Alpha_R2;
end

% % Load the data into one matrix
SampleSizeHC = 17;
SampleSizeDZ = 17;
NEEGPoints1 = 130; % ERT task
NEEGPoints2 = 110; % GRT task - positive reward
NEEGPoints3 = 70;  % GRT task - positive feedback
CnctDim = 34;
NTests = 3;
NPointsHC = NEEGPoints1*SampleSizeHC;
NPointsDZ = NEEGPoints1*SampleSizeDZ;
TotalNPoints = (SampleSizeHC+SampleSizeDZ)*(NTests*NEEGPoints1+NEEGPoints2+NEEGPoints3);
% control group
DyMatAll = zeros(CnctDim, CnctDim, TotalNPoints);
for subjId = 1: SampleSizeHC
    DyMatAll(:,:, (1+NEEGPoints1*(subjId-1)):(NEEGPoints1*(subjId))) = abs(ControlA_N{1,subjId});
    DyMatAll(:,:, (1+NEEGPoints1*(subjId-1) + NPointsHC):(NEEGPoints1*(subjId)) + NPointsHC) = abs(ControlA_M{1,subjId});
    DyMatAll(:,:, (1+NEEGPoints1*(subjId-1) + 2*NPointsHC):(NEEGPoints1*(subjId) + 2*NPointsHC)) = abs(ControlA_R{1,subjId});
end

% Disease group

for subjId = 1: SampleSizeDZ
    DyMatAll(:,:, (1+NEEGPoints1*(subjId-1) + 3*NPointsDZ):(NEEGPoints1*(subjId) + 3*NPointsDZ)) = abs(DiseaseA_N{1,subjId});
    DyMatAll(:,:, (1+NEEGPoints1*(subjId-1) + 4*NPointsDZ):(NEEGPoints1*(subjId) + 4*NPointsDZ)) = abs(DiseaseA_M{1,subjId});
    DyMatAll(:,:, (1+NEEGPoints1*(subjId-1) + 5*NPointsDZ):(NEEGPoints1*(subjId) + 5*NPointsDZ)) = abs(DiseaseA_R{1,subjId});
end

%  award
NEEGPoints2 = 70;
NPointsHC2 = NEEGPoints2*SampleSizeHC;


for subjId = 1: SampleSizeDZ
    DyMatAll(:,:, (1+NEEGPoints2*(subjId-1) + 6*NPointsDZ):(NEEGPoints2*(subjId) + 6*NPointsDZ)) = abs(ControlA_CW{1,subjId});
    DyMatAll(:,:, (1+NEEGPoints2*(subjId-1) + 6*NPointsDZ+NPointsHC2):(NEEGPoints2*(subjId) + 6*NPointsDZ+NPointsHC2)) = abs(DiseaseA_CW{1,subjId});
end

%feedback
NEEGPoints3 = 110;
NPointsHC3 = NEEGPoints3*SampleSizeHC;


for subjId = 1: SampleSizeDZ
    DyMatAll(:,:, (1+NEEGPoints3*(subjId-1) + 6*NPointsDZ+2*NPointsHC2):(NEEGPoints3*(subjId) + 6*NPointsDZ+2*NPointsHC2)) = abs(ControlA_AR{1,subjId});
    DyMatAll(:,:, (1+NEEGPoints3*(subjId-1) + 6*NPointsDZ+2*NPointsHC2+NPointsHC3):(NEEGPoints3*(subjId) + 6*NPointsDZ+2*NPointsHC2+NPointsHC3)) = abs(DiseaseA_AR{1,subjId});
end
% %clear diagnoal
% for i=1:TotalNPoints
%   DyMatAll(i,i)=0;
% end
%% Run dissimilarity here
% compute Euclidean distance from each pair of connectomes
X = reshape(DyMatAll,[CnctDim*CnctDim,TotalNPoints]);

% L2
D = sum(X .^ 2);
DyDistAll = real(sqrt(bsxfun(@plus, D.', D) - (2 * (X.' * X))));
DyDistAll(1:length(DyDistAll)+1:end) = 0;
clear X D % clear memory from large arrays
%% Run NDR here
Edim=10;
%
n=30;
[IsomapXYZ, dumpAll]=compute_mapping(DyDistAll ,'Isomap', Edim, n);

% load('conn_comp.mat');
%%
IsomapAll=zeros(length(dumpAll.conn_comp),3);
IsomapAll(dumpAll.conn_comp,:)=IsomapXYZ(:,1:3);
%%
% %% plot output
% figure;
% label=zeros(3,2600);
% for k=1:3
%     label(k,:)=1+2600*(k-1):2600*k;
% end
% m_size=2;
%
% marker_hue=['g','m','b','g','m','b',];
% marker_style=['o','o','o','o','o','o'];
% for k=1:3
%     plot3(IsomapAll(label(k,:),1), IsomapAll(label(k,:),2), IsomapAll(label(k,:),3), marker_style(k),'MarkerEdgeColor',marker_hue(k),'MarkerFaceColor',marker_hue(k),'MarkerSize',m_size);
%
%     grid on;
%
%     axis equal;
%     hold on
% end
%
% %%  plot average track using the out of sample
% % mainHC=21:40;
% % reappHC=41:60;
% % neuDZ=61:80;
% % mainDZ=81:100;
% % reappDZ=101:120;
%
% %AllDist = reshape(DyDistAll,[120,NEEGPoints,TotalNPoints]);
% %MeanTrack=zeros(2*NTests,NEEGPoints,TotalNPoints);
% % for i=6%:(2*NTests)
% %     MeanTrack(i,:,:) = squeeze(mean(AllDist(1+(i-1)*20:i*20,:,:),1));
% % end
% % tic;
% % %OutOfSample_Mean = zeros(2*NTests,NEEGPoints,Edim);
% % %  parpool(4);
% % for i = 6%:(2*NTests)
% %     OutOfSample_Mean(i,:,:) = out_of_sample( squeeze(MeanTrack(i,:,:)), dumpAll);
% % end
% %  %parpool('close');
% % toc
% %%
% % line_style=['-o','-o','-o','-*','-*','-*',];
% % line_size=4;
% % figure;
% % for i=1:6
% %     subplot(2,3,i)
% %     plot3( OutOfSample_Mean(i,: ,1), OutOfSample_Mean(i,: ,2), OutOfSample_Mean(i,: ,3), line_style(i),'MarkerEdgeColor',marker_hue(i),'MarkerFaceColor',marker_hue(i),'MarkerSize',line_size);
% %     axis equal;
% %     grid on;
% %     hold on;
% % end
% % hold off
%
%

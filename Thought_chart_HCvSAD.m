%% This is the main processing code for dissimilarity based Thought Chart
% for subject conducting ERT task, there are two groups HC vs SAD subjects
% in three conditions: neutral, maintain and reappraise. Input are dynmaic
% connectivity matrix each in size 34*34*50*130, and foi averaged matrix are
% taken for dissimilarity matrix construction.
% Noted that Neutral is alway marked in green
% maintain is in magenta and reappraise is in blue
% script write in 05/29/2016 by Mengqi Xing

clear all
%% loading and finding the foi average for 2 groups
ConnectC={'1','6','9','13','15','17','19','22','35','39','42','43','44','47','49','51','52','53','54','57'};
ConnectP={'1','2','6','8','17','21','23','30','37','40','50','118','124','127','138','140','143','149','153','167'};
foi=5;%define the frequency of interest

%control
for i=1:length(ConnectC)
    
    
    RestC=['WPLI_Reappraise_HC' ConnectC{i} '.mat'];
    R=cell2mat(struct2cell(load(RestC)));
    Alpha_R(:,:,:)=abs(mean(R(:,:,foi,1:130),3));
    ControlA_R{i}=Alpha_R;
    
    RestC=['WPLI_Neutral_HC' ConnectC{i} '.mat'];
    N=cell2mat(struct2cell(load(RestC)));
    Alpha_N(:,:,:)=abs(mean(N(:,:,foi,1:130),3));
    ControlA_N{i}=Alpha_N;
    
    RestC=['WPLI_Maintain_HC' ConnectC{i} '.mat'];
    M=cell2mat(struct2cell(load(RestC)));
    Alpha_M(:,:,:)=abs(mean(M(:,:,foi,1:130),3));
    ControlA_M{i}=Alpha_M;
end

% disease
for i=1:length(ConnectP)
    
    
    RestP=['WPLI_Reappraise_DZ' ConnectP{i} '.mat'];
    R=cell2mat(struct2cell(load(RestP)));
    Alpha_R(:,:,:)=mean(R(:,:,foi,1:130),3);
    DiseaseA_R{i}=Alpha_R;
    
    RestP=['WPLI_Neutral_DZ' ConnectP{i} '.mat'];
    N=cell2mat(struct2cell(load(RestP)));
    Alpha_N(:,:,:)=mean(N(:,:,foi,1:130),3);
    DiseaseA_N{i}=Alpha_N;
    
    RestP=['WPLI_Maintain_DZ' ConnectP{i} '.mat'];
    M=cell2mat(struct2cell(load(RestP)));
    Alpha_M(:,:,:)=mean(M(:,:,foi,1:130),3);
    DiseaseA_M{i}=Alpha_M;
end


%% Load the data into one matrix
SampleSizeHC= 20  ;
SampleSizeDZ= 20  ;
%control group


for NSubj = 1: SampleSizeHC
    temp = ControlA_N{1,NSubj};
    for i=1:130
        DyMatAll(:,:, i+130*(NSubj-1)) = abs( temp(:, :,  i) ) ;
    end
    
    temp = ControlA_M{1,NSubj};
    for i=1:130
        DyMatAll(:,:, i+130*(NSubj-1) +2600) = abs( temp(:, :, i) ) ;
    end
    
    temp = ControlA_R{1,NSubj};
    for i=1:130
        DyMatAll(:,:, i+130*(NSubj-1)+5200) = abs( temp(:, :,  i) ) ;
    end
end

% Disease group

for NSubj = 1: SampleSizeDZ
    temp = DiseaseA_N{1,NSubj};
    for i=1:130
        DyMatAll(:,:, i+130*(NSubj-1)+7800 ) = abs( temp(:, :,  i) ) ;
    end
    
    temp = DiseaseA_M{1,NSubj};
    for i=1:130
        DyMatAll(:,:, i+130*(NSubj-1)+10400 ) = abs( temp(:, :,  i) ) ;
    end
    
    temp = DiseaseA_R{1,NSubj};
    for i=1:130
        DyMatAll(:,:, i+130*(NSubj-1)+13000) = abs( temp(:, :,  i) ) ;
    end
end

%
% clear diagnoal

for i=1:40*3*130
    for j=1:34
        DyMatAll(j,j,i)=0;
    end
end
dim=40*3*130;
%% Run dissiilarity here
DyDistAll=zeros(dim, dim);

for i=1:dim
    i
    for j=(i+1): dim
        diff = DyMatAll(:,:, i) -DyMatAll(:,:,j);
        DyDistAll(i, j) = sqrt( sum( sum( diff.*diff ) ) );
        DyDistAll(j, i) = DyDistAll(i,j);
    end
end
%% Run NDR here
Edim=10;

n=30;

[IsomapXYZ dumpAll]=compute_mapping(DyDistAll ,'Isomap', Edim,n);
load('conn_comp.mat');
%
IsomapAll=zeros(length(conn_comp),3);
IsomapAll(conn_comp,:)=IsomapXYZ(:,2:4);
m_size=2;

%% plot output
figure;
label=zeros(6,2600);
for k=1:6
    label(k,:)=1+2600*(k-1):2600*k;
end
m_size=2;

marker_hue=['g','m','b','g','m','b',];
marker_style=['o','o','o','*','*','*',];
for k=1:6
    plot3(IsomapAll(label(k,:),1), IsomapAll(label(k,:),2), IsomapAll(label(k,:),3), marker_style(k),'MarkerEdgeColor',marker_hue(k),'MarkerFaceColor',marker_hue(k),'MarkerSize',m_size);
    
    grid on;
    
    axis equal;
    hold on
end

%%  plot average track using the out of sample 
mainHC=21:40;
reappHC=41:60;
neuDZ=61:80;
mainDZ=81:100;
reappDZ=101:120;

AllDist=reshape(DyDistAll,[120,130,15600]);
MeanTrack=zeros(6,15600);
for i=1:6
MeanTrack(i,:)=mean(AllDist(1+(i-1)*20:i*20,:,:),1);
end

OutOfSample_Mean= out_of_sample( MeanTrack, dumpAll);
line_style=['-o','-o','-o','-*','-*','-*',];
line_size=4;
for i=1:6
plot3( OutOfSample_Mean(: ,1), OutOfSample_Mean(: ,2), OutOfSample_Mean(: ,3), line_style(i),'MarkerEdgeColor',marker_hue(i),'MarkerFaceColor',marker_hue(i),'MarkerSize',line_size);
axis equal;
grid on;
hold on;
end
hold off



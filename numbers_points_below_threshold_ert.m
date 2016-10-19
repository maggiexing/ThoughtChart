%% This is the code used for calculating the CDF for the Thought Chart

%% change the threshold here
threshold=-1900:50:300;

N_HC=IsomapAll(1:2600,1);
M_HC=IsomapAll(2601:5200,1);
R_HC=IsomapAll(5201:7800,1);
N_DZ=IsomapAll(7801:10400,1);
M_DZ=IsomapAll(10401:13000,1);
R_DZ=IsomapAll(13001:15600,1);


% healthy
for i=1:length(threshold)
    count_NHC=0;
for j=1:2600
        if N_HC(j)<(threshold(i)) 
        count_NHC=count_NHC+1;
    end
end
Total_NHC(i)=count_NHC/2600;
end

for i=1:length(threshold)
    count_MHC=0;
for j=1:2600
        if M_HC(j)<(threshold(i))
        count_MHC=count_MHC+1;
    end
end
Total_MHC(i)=count_MHC/2600;
end

for i=1:length(threshold)
    count_RHC=0;
for j=1:2600
        if R_HC(j)<(threshold(i)) 
        count_RHC=count_RHC+1;
    end
end
Total_RHC(i)=count_RHC/2600;
end

% disease
for i=1:length(threshold)
    count_NDZ=0;
for j=1:2600
        if N_DZ(j)<(threshold(i)) 
        count_NDZ=count_NDZ+1;
    end
end
Total_NDZ(i)=count_NDZ/2600;
end

for i=1:length(threshold)
    count_MDZ=0;
for j=1:2600
        if M_DZ(j)<(threshold(i))
        count_MDZ=count_MDZ+1;
    end
end
Total_MDZ(i)=count_MDZ/2600;
end

for i=1:length(threshold)
    count_RDZ=0;
for j=1:2600
        if R_DZ(j)<(threshold(i)) 
        count_RDZ=count_RDZ+1;
    end
end
Total_RDZ(i)=count_RDZ/2600;
end
figure();
plot(-threshold,-Total_NHC,'-g')
hold on
plot(-threshold,-Total_MHC,'-m')
hold on
plot(-threshold,-Total_RHC,'-b')
hold on
plot(-threshold,-Total_NDZ,'--g')
hold on
plot(-threshold,-Total_MDZ,'--m')
hold on
plot(-threshold,-Total_RDZ,'--b')
xlabel('Threshold')
ylabel('Numbers of the points below the threshold')
legend('Neutral HC','Maintain HC','Reappraise HC','Neutral DZ','Maintain DZ','Reappraise DZ')
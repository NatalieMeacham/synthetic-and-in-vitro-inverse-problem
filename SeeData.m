%Visualize all replicates of a specific dosage for a specific dataset.
%Note that both the dataset choice and the relevant figure legend need to be
%uncommented.

load('MONOCLONAL_DATA.mat');
load('Mixture_Data.mat')

%Uncomment T and Tstring for desired dataset:

% T=RESISTANT_250_BF; 
%     Tstring="R250";
% T=RESISTANT_500_BF; 
%     Tstring="R500";
% T=SENSITIVE_500_BF; 
%     Tstring="S500";
%T=SENSITIVE_1000_BF; 
    %Tstring="S100";
%T=BF_11;
    %Tstring="BF11";
% T=BF_12;
%     Tstring="BF12";
T=BF_21;
    Tstring="BF21";
%T=BF_41;
    %Tstring="BF41";

%choose concentration
chooseconc=5;
concstr=string(chooseconc);

sz = size(T);
tsize=sz(3);
concsize=sz(2);
repsize=sz(1);
tvec=linspace(1,tsize,tsize);
concvec=linspace(1,concsize,concsize);
meanvec = zeros(tsize,1);
meanmat = zeros(concsize,tsize);
concnum=length(concvec);
tspan=linspace(9,48,14);

for i=tvec
    for j=concvec
        m=mean(T(:,j,i),'omitnan');
        meanmat(j,i)=m;
    end
end

%cross-section of the data for a given concentration
matrix1=zeros(sz(1),sz(3));
for r=1:sz(1)
    for t=1:sz(3)
        matrix1(r,t)=T(r, chooseconc,t );
    end
end

figure
for r=1:sz(1)
    plot(tspan,matrix1(r,:),'o','LineWidth',2)
    hold on
end

plot(tspan,meanmat(chooseconc,:),'LineWidth',2)
%monoclonal data legend
%legend('Replicate 1','Replicate 2','Replicate 3','Replicate 4','Replicate 5','Replicate 6','Replicate 7','Mean','Location','Northwest')
%mixture data legend
legend('Replicate 1','Replicate 2','Replicate 3','Replicate 4','Replicate 5','Replicate 6','Replicate 7','Replicate 8','Replicate 9','Replicate 10','Replicate 11','Replicate 12','Replicate 13','Replicate 14','Mean','Location','Northwest','FontSize',14)
xlabel('Time')
ylabel('Number of Cells')
title(strcat('Dataset',{' '},Tstring,{' '},'for Dosage ',{' '},concstr))
set(gca,"FontSize",20)


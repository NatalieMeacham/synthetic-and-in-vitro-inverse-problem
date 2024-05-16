load('MONOCLONAL_DATA.mat');
load('Mixture_Data.mat')
%T=RESISTANT_250_BF; %first row has NaNs, use concvec=linspace(2,concsize,concsize-1);
T=RESISTANT_500_BF; %last row has NaN, use concvec=linspace(1,concsize-1,concsize-1);
%T=SENSITIVE_500_BF; %no NaN
%T=SENSITIVE_1000_BF; %last row has NaN, use concvec=linspace(1,concsize-1,concsize-1);
T=BF_11;
sz = size(T);
tsize=sz(3);
concsize=sz(2);
repsize=sz(1);
tvec=linspace(1,tsize,tsize);
concvec=linspace(1,concsize,concsize);
%concvec=linspace(2,concsize,concsize-1);
%concvec=linspace(1,concsize-1,concsize-1)
meanvec = zeros(tsize,1);
meanmat = zeros(concsize,tsize);
concnum=length(concvec);
tspan=linspace(9,48,14);

for i=tvec
    for j=concvec
        m=mean(T(:,j,i),'omitnan');
    %m = mean(T(i,:),'omitnan');
    %meanvec(i)=m;
        meanmat(j,i)=m;
    end
end

%do for loop to assemble vector of data that depicts each replicate over
%time (actually a matrix)
chooseconc=5;

%matrix for a given concentration with rows being each replicate and cols
%being each time step
%so this matrix is a cross-section of the data for a given concentration
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
%plot(meanmat(chooseconc,:),'LineWidth',2)
plot(tspan,meanmat(chooseconc,:),'LineWidth',2)
legend('Replicate 1','Replicate 2','Replicate 3','Replicate 4','Replicate 5','Replicate 6','Replicate 7','Mean','Location','Northwest')
legend('Replicate 1','Replicate 2','Replicate 3','Replicate 4','Replicate 5','Replicate 6','Replicate 7','Replicate 8','Replicate 9','Replicate 10','Replicate 11','Replicate 12','Replicate 13','Replicate 14','Mean','Location','Northwest','FontSize',8)
xlabel('Time')
ylabel('Number of Cells')
set(gca,"FontSize",20)

%see data for a given concentration at all time points and replicates
matrix1
meanmat(chooseconc,:)

% figure
% dosevec= [0, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];
% plot(1:11,dosevec,'d','LineWidth',2)
% %title('$$Q \geq \frac{I_h H}{I_h H+I_z C}, b_1 \geq b_2$$','interpreter','latex')
% xlabel('Concentration Number')
% string=strcat('Dosage Level in','\mu M','interpreter','latex')
% string=strcat('Dosage Level in ', {' '},'\muM')
% %ylabel('Dosage Level')
% ylabel(string)
% set(gca,"FontSize",20)
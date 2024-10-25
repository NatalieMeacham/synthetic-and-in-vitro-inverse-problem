clc
clear all
%The goal of this script is to update GLSInverseScriptData by using the
%drug concentration values given by PhenoPop to individualize the growth
%and death rates for every concentration
%Note that there are only 10 concentration doses given, so I am assuming
%that the 11th is a control/0 dose. 

format long

choosedata=['BF12']; %choose between R250, R500, S500, and S1000
%or BF11,BF12,BF21,BF41
%rho=0.3; %originially chosen from phenopop to be 0.3
%k=0.15; %originially chosen arbitrarily to be 0.15
%rho=0.0633;
%rho=0.0634
%k=0.047;
%kmax=0.0047;
kmax=0.0077; %regular
 %kmax = 0.0019; %when we have death term incl rho?
kmin=0;
%k=0;
%rho=0.0653;
%rho=0.0682;
rhomax=0.0682;
%rhomax=0.0692
%rhomax=0.0721;
%rhomin=0.0237;
rhomin=0.0375;
%rhomin=0.0241
%rhomin=0.08;
%rhomin=0;

if choosedata=='R250'
    rhomin=0.0241;
    rhomax=0.0692; %from first row
    kmin=0;
    kmax=0.0019; %just chose lower option
elseif choosedata=='R500'
    rhomin=0.0375;
    rhomax=0.0682; 
    kmin=0;
    kmax=0.0019; %just chose lower option
elseif choosedata=='S500'
    rhomin=0.0241; %just chose lower option
    rhomax=0.0692; %just chose higher option 
    kmin=0;
    kmax=0.0019; %from 10th row
elseif choosedata=='S100'
    rhomin=0.0241; %just chose lower option
    rhomax=0.0692; %just chose higher option 
    kmin=0;
    kmax=0.0077; %from 9th row
    %kmax=0.0047 %if we have rho in death term
else
    disp('??')
end
dosevec= [0, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];
dosevecshort= [0.00, 0.03,0.05, 0.12, 0.22, 0.34,0.46,1.19,2.44,3.71,5.01];
dosevecstring= ["0.00", "0.03","0.05", "0.12", "0.22", "0.34","0.46","1.19","2.44","3.71","5.01"];
%rhovec=rho.*(1./(1 + dosevec));
%rhovec=rho*ones(length(dosevec),1);
%rhovec=rhomax*ones(length(dosevec),1);

%rhovec=zeros(length(dosevec),1);
%rhovec=rho.*(dosevec./dosevec(11));
%rhovec=rho*ones(length(dosevec),1) - rho.*(dosevec./dosevec(11))';
%kvec=k.*0.5*dosevec; %now apply these values by concentration in each for loop
%kvec=k.*dosevec;
kvec=kmax.*(dosevec./dosevec(11));
%kvec=linspace(0,k,11);
%kvec=k.*(ones(length(dosevec)) + dosevec);
%kvec=k.*(1 + (dosevec./dosevec(11)));
%kvec=k.*ones(length(dosevec),1) + (dosevec./dosevec(11));

% %scale rhovec from dosevec:
% rhodiff=rhomax-rhomin;
% scaledosevec=zeros(length(dosevec),1);
% rhovec=zeros(length(dosevec),1);
% for i=1:length(dosevec)
%     scaledosevec(i)=dosevec(i)/dosevec(11)*(rhodiff) + rhomin;
% end
% 
% %the larger concentrations should have higher dosages->lower growth rate
% for i=1:length(dosevec)
%     rhovec(i)=scaledosevec(12-i);
% end

%scale rhovec from dosevec 
rhodiff=rhomax-rhomin;
scaledosevec=zeros(length(dosevec),1);
rhovec=zeros(length(dosevec),1);
for i=1:length(dosevec)
    scaledosevec(i)=1 - dosevec(i)/dosevec(11);
end

%the larger concentrations should have higher dosages->lower growth rate
for i=1:length(dosevec)
    rhovec(i)=scaledosevec(i)*rhodiff + rhomin;
end

tinit=9;
tfinal=48;
tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

[meanmatnorm,concvec,errmat]=ScaleData(choosedata);

[meanmatnorm,isnanmat, tspanvec, tspanmat,concvecA]=DropNaNsFn(concvec,meanmatnorm,tspan);

%find AIC vals for each concentration for each number of rpoints
%rpointsvec=6:2:24;
rpointsvec=4:1:30;
AICmat=zeros(length(concvecA),length(rpointsvec));
for a=1:length(concvecA)
    data=meanmatnorm(a,1:tspanvec(a));
    %keyboard
    y0=data(1);
    rho=rhovec(concvecA(a));
    k=kvec(concvecA(a));
    %keyboard
    tspanS=tspanmat(a,1:tspanvec(a));
    for i=1:length(rpointsvec)
    [~,~,AIC_GLS,~,~]=GLSInverseFnNData(rho,k,y0,tspanS,rpointsvec(i),data);
    AICmat(a,i)=AIC_GLS;
    end
end

%create a matrix that fills in the optweight corresponding to the lowest
%AIC value for each conc
%modified with dropped nans
%concnum=length(concvec);
concnum=length(concvecA);
AICminvec=zeros(concnum,1);
AICcheck=zeros(concnum,1);
optmat=zeros(concnum,length(rpointsvec));
for a=1:length(concvecA)
    data=meanmatnorm(a,1:tspanvec(a))
    %keyboard
    tspanS=tspanmat(a,1:tspanvec(a));
    %keyboard
    rho=rhovec(concvecA(a));
    k=kvec(concvecA(a));
    [M,I]=min(AICmat(a,:));
    AICminvec(a)=rpointsvec(I);
    %data=meanmatnorm(a,:);
    %y0=meanmatnorm(a,1);
    y0=data(1);
    %tspan=tspan(1:tspanvec(a));
    [optweight,~,~,~,~]=GLSInverseFnNData(rho,k,y0,tspanS,AICminvec(a),data);
    optmat(a,1:AICminvec(a))=optweight;
end

concvecS=concvecA;

%create legend for plots
for a=1:length(concvecS)
    for iter=1:length(concvecS)
        Legend{iter}=strcat('Conc.',' ', num2str(concvecS(iter)));
        Legend2{iter}=strcat('d=',num2str(dosevecshort(concvecS(iter))),'  \muM');
        Legend3{iter}=strcat('{\it d}=',dosevecstring(concvecS(iter)));
    end
end

%clr = cool(length(concvecS));
clr = hsv(length(concvecS));

%create matrices of optimal grids and optimal dists 
trymatoptrsgrid=zeros(length(concvecS),max(AICminvec));
for a=1:length(concvecS)
    optrpoints=AICminvec(a);
    optrsgrid=linspace(0,1,optrpoints);
    trymatoptrsgrid(a,1:optrpoints)=optrsgrid;
end

%create matrix and compute cdfs for all dists
cdfRmat=zeros(length(concvecS),max(AICminvec));
for a=1:length(concvecS)
    optrpoints=AICminvec(a);
    optrsgrid=linspace(0,1,optrpoints);
    partialsumR=0;
    for i=1:optrpoints
        partialsumR=sum(optmat(a,1:i));
        cdfRmat(a,i)=partialsumR;
    end
end

%concvecS=concvec;
trymatoptrsgridS=trymatoptrsgrid;
optmatS=optmat; 
AICminvecS=AICminvec;
cdfRmatS=cdfRmat; 
meanmatnormS=meanmatnorm;

%PLOT PDFS CORRECTLY
figure
for a=1:length(concvecS)
    hold on
    %plot(trymatoptrsgridS(a,1:AICminvecS(a)),optmatS(a,1:AICminvecS(a)),'LineWidth',2,'Color',clr(a,:))
    %plot(trymatoptrsgridS(a,1:AICminvecS(a)),optmatS(a,1:AICminvecS(a)),'--o','MarkerSize',8,'LineWidth',2,'Color',clr(a,:))
    %plot(trymatoptrsgridS(a,1:AICminvecS(a)),optmatS(a,1:AICminvecS(a)),'--o','MarkerSize',8,'LineWidth',2,'MarkerEdgeColor',clr(a,:),'MarkerFaceColor',clr(a,:))
    plot(trymatoptrsgridS(a,1:AICminvecS(a)),optmatS(a,1:AICminvecS(a)),'--o','Color',clr(a,:),'MarkerFaceColor',clr(a,:),'MarkerSize',8,'LineWidth',2)
    
end
xlabel('Sensitivity to Treatment {\it s}')
%ylabel('Recovered Prop. of Pop.') %original version 
ylabel('Recovered Proportion of Population')
legend(Legend, 'FontSize',12)
legend(Legend2, 'FontSize',12)
legend(Legend3, 'FontSize',12,'Location','northwest')
set(gca,"FontSize",20)

%plot cdfs again w stair fn
figure
for a=1:length(concvecS)
    hold on
    %X = linspace(0,4*pi,40);
    clr = hsv(length(concvecS));
    X=trymatoptrsgridS(a,1:AICminvecS(a));
    %Y = sin(X);
    Y=cdfRmatS(a,1:AICminvecS(a));
    stairs(X,Y,'--o','Color',clr(a,:),'MarkerFaceColor',clr(a,:),'MarkerSize',8,'LineWidth',2)
    ylim([0 1])
    xlabel('Sensitivity to Treatment {\it s}')
%ylabel('Cumulative Recovered Prop. of Pop.') %original version
ylabel('Cumulative Recovered Proportion')
legend(Legend3, 'FontSize',12,'Location','northwest')
set(gca,"FontSize",20)
end

%plot data vs weightedsol for non-nan concentrations NEW VERSION
wsolmat=zeros(length(concvecS),length(tspan));
figure
for a=1:length(concvecS)
    %define rsgrid for each a
    rho=rhovec(concvecA(a));
    k=kvec(concvecA(a));
    %data=meanmatnormS(a,:);
    data=meanmatnorm(a,1:tspanvec(a));
    %keyboard
    tspanS=tspanmat(a,1:tspanvec(a));
    %keyboard
    y0=data(1);
    optrpoints=AICminvecS(a);
    optrsgrid=linspace(0,1,optrpoints);
    optweight=optmatS(a,1:AICminvecS(a));
   % keyboard
    [~, ~,weightedsol] = ForwardFunctionN(optrsgrid, optweight, rho, k, y0, tspanS);
    wsolmat(a,1:length(tspanS))=weightedsol;
    clr = hsv(length(concvecS));    %# LINES colormap
    %plot(tspanS,data,'d','Color',clr(a,:),'MarkerFaceColor',clr(a,:),'MarkerSize',8,'LineWidth',2,'HandleVisibility','off') %version without error bar
    errorbar(tspanS,data,errmat(a,1:tspanvec(a)),'.','MarkerSize',8,'LineWidth',2,'Color',clr(a,:),'HandleVisibility','off')
    hold on
    plot(tspanS,weightedsol,'LineWidth',2,'Color',clr(a,:))
    hold on
    % xlabel('Time (Hours)')
    % ylabel('Tumor Growth Data and Model Fit')
    % legend(Legend3, 'FontSize',12,'Location','northwest')
    % set(gca,"FontSize",20)
end
xlabel('Time (Hours)')
    ylabel('Tumor Growth Data and Model Fit')
    legend(Legend3, 'FontSize',12,'Location','northwest')
    set(gca,"FontSize",20)

%compute expected value of the recovered dist
evvec=zeros(length(concvecS),1);
for a=1:length(concvecS)
    optrpoints=AICminvecS(a);
    optrsgrid=linspace(0,1,optrpoints);
    optweight=optmatS(a,1:AICminvecS(a));
    evvec(a)=sum(optrsgrid.*optweight);
end
% %%
% %plot data with error bars
% %figure
% for a=1:length(concvecS)
%     figure
%     errorbar(tspanmat(a,1:tspanvec(a)),meanmatnorm(a,:),errmat(a,:),'Color',clr(a,:))
%     title(num2str(a))    %hold on
% end 
% figure
% plot(concvecS,evvec,'o','LineWidth',2)
% hold on
% xlabel('Concentration Number')
% ylabel('Mean Sensitivity')
% set(gca,"FontSize",20)
% 
% figure
% plot(dosevec,evvec,'o','LineWidth',2)
% hold on
% xlabel('Dosage Values')
% ylabel('Mean Sensitivity')
% set(gca,"FontSize",20)

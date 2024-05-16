clc
clear all
%The goal of this script is to update GLSInverseScriptData by using the
%drug concentration values given by PhenoPop to individualize the growth
%and death rates for every concentration
%Note that there are only 10 concentration doses given, so I am assuming
%that the 11th is a control/0 dose. 

format long

choosedata=['S500']; %choose between R250, R500, S500, and S1000
%or BF11,BF12,BF21,BF41
%rho=0.3; %originially chosen from phenopop to be 0.3
%k=0.15; %originially chosen arbitrarily to be 0.15
%rho=0.0633;
%rho=0.0634
%k=0.047;
%kmax=0.0047;
kmax=0.0077;
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
    % %v2
    % rhomin=0;
    % rhomax=0.0452; %this made it worse lol
    kmin=0;
    kmax=0.0019; %from 10th row
elseif choosedata=='S100'
    rhomin=0.0241; %just chose lower option
    rhomax=0.0692; %just chose higher option 
    kmin=0;
    kmax=0.0077; %from 9th row
else
    disp('??')
end
dosevec= [0, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];
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

%scale rhovec from dosevec BETTER SCALE?:
rhodiff=rhomax-rhomin;
scaledosevec=zeros(length(dosevec),1);
rhovec=zeros(length(dosevec),1);
for i=1:length(dosevec)
    %scaledosevec(i)=dosevec(i)/dosevec(11);
    scaledosevec(i)=1 - dosevec(i)/dosevec(11);
end

%the larger concentrations should have higher dosages->lower growth rate
for i=1:length(dosevec)
    %rhovec(i)=scaledosevec(12-i)*rhodiff + rhomin;
    rhovec(i)=scaledosevec(i)*rhodiff + rhomin;
end

% kdiff=kmax-kmin;
% scaledoseveck=zeros(length(dosevec),1);
% kvec=zeros(length(dosevec),1);
% for i=1:length(dosevec)
%     %scaledosevec(i)=dosevec(i)/dosevec(11);
%     scaledoseveck(i)=1 - dosevec(i)/dosevec(11);
% end
% %scaledoseveck=flip(scaledoseveck);
% 
% %the larger concentrations should have higher dosages->lower growth rate
% for i=1:length(dosevec)
%     %rhovec(i)=scaledosevec(12-i)*rhodiff + rhomin;
%     kvec(i)=scaledoseveck(i)*kdiff + kmin;
% 
% end
% kvec=flip(kvec);
% 
% kvec=linspace(kmin,kmax,11);
%try more linear rhovec:
%rhovec=linspace(rhomax,rhomin,11);
% 
% keyboard

%import data 
% load('MONOCLONAL_DATA.mat');
tinit=9;
tfinal=48;
tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

[meanmatnorm,concvec]=ScaleData(choosedata);
%keyboard
%meanmatnorm
%keyboard

%OLD MAXVAL (depends on the max of the data matrix)
%maxval = 10*max(meanmat,[],'all');

% %normalize the data
% meanmatnorm=zeros(size(meanmat));
% for a=concvec
%     for t=tvec
%         s=meanmat(a,t)/maxval; %multiply by 10 because data still growing?
%         meanmatnorm(a,t)=s;
%     end 
% end

%keyboard

% %replace any NaNs row value with 1s
% TFvec=zeros(length(concvec),1);
% for a=concvec
%     TFvec(a)=anynan(meanmatnorm(a,:))
%     %TFvec=[TFvec, TF]
% end
% 
% for a=concvec
%     if TFvec(a)==1
%         meanmatnorm(a,:)=1;
%         disp('updated meanmatnorm with ones')
%     else 
%         disp('chill')
%     end
% end
%keyboard
%[meanmatnorm,tspanvec]=DropNaNsFn(concvec,meanmatnorm,tspan);
%[meanmatnorm,isnanmat, tspanvec, tspanmat]=DropNaNsFn(concvec,meanmatnorm,tspan)
[meanmatnorm,isnanmat, tspanvec, tspanmat,concvecA]=DropNaNsFn(concvec,meanmatnorm,tspan)
%keyboard
%keyboard
%keyboard 
%tspanvec;
%meanmatnorm
%keyboard

% tinit=9;
% tfinal=48;
% tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

%OLD VERSION BEFORE DROPPING NAN ROWS
% %find AIC vals for each concentration for each number of rpoints
% rpointsvec=6:2:24;
% AICmat=zeros(length(concvec),length(rpointsvec));
% for a=concvec
%     data=meanmatnorm(a,1:tspanvec(a))
%     %keyboard
%     y0=meanmatnorm(a,1);
%     rho=rhovec(a);
%     k=kvec(a);
%     tspanS=tspanmat(a,1:tspanvec(a));
%     for i=1:length(rpointsvec)
%     [~,~,AIC_GLS,~,~]=GLSInverseFnNData(rho,k,y0,tspanS,rpointsvec(i),data);
%     AICmat(a,i)=AIC_GLS;
%     end
% end

%find AIC vals for each concentration for each number of rpoints
rpointsvec=6:2:24;
AICmat=zeros(length(concvecA),length(rpointsvec));
for a=1:length(concvecA)
    data=meanmatnorm(a,1:tspanvec(a))
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

%TFvec=[0 0 0 0 0 0 0 0 0 0 0];
%keyboard

% %define concvecS as concvec without any nan concentrations
% for a=concvec
%     %define rsgrid for each a
%     optrpoints=AICminvec(a);
%     optrsgrid=linspace(0,1,optrpoints);
%     clr = cool(length(concvec));
%     if sum(TFvec)==1
%         if TFvec(a)==1
%             concvecS=concvec;
%             concvecS(a)=[];
%         else
%             %disp('fsdfsdf')
%         end
%     else
%         concvecS=concvec;
%     end
% end

concvecS=concvecA;

%create legend for plots
for a=1:length(concvecS)
    for iter=1:length(concvecS)
        Legend{iter}=strcat('Conc.',' ', num2str(concvecS(iter)));
    end
end

clr = cool(length(concvecS));

% %create matrices of optimal grids and optimal dists 
% trymatoptrsgrid=zeros(concnum,max(AICminvec));
% for a=concvec
%     optrpoints=AICminvec(a);
%     optrsgrid=linspace(0,1,optrpoints);
%     trymatoptrsgrid(a,1:optrpoints)=optrsgrid;
% end

%create matrices of optimal grids and optimal dists 
trymatoptrsgrid=zeros(length(concvecS),max(AICminvec));
for a=1:length(concvecS)
    optrpoints=AICminvec(a);
    optrsgrid=linspace(0,1,optrpoints);
    trymatoptrsgrid(a,1:optrpoints)=optrsgrid;
end

% %create matrix and compute cdfs for all dists
% cdfRmat=zeros(length(concvec),max(AICminvec));
% for a=concvec
%     optrpoints=AICminvec(a);
%     optrsgrid=linspace(0,1,optrpoints);
%     partialsumR=0;
%     for i=1:optrpoints
%         partialsumR=sum(optmat(a,1:i));
%         cdfRmat(a,i)=partialsumR;
%     end
% end

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

% %cut of NaN row for all relevant matrices and vectors if necessary
% for a=concvec
%     if sum(TFvec)==1
%         if TFvec(a)==1
%             %A=strcat('we are not plotting concentration #',num2str(a), ' due to NaNs');
%             %disp('we are not plotting %i', a)
%             %disp(A)
%             % concvecS=concvec;
%             % concvecS(a)=[];
%             trymatoptrsgridS=trymatoptrsgrid;
%             trymatoptrsgridS(a,:)=[];
%             optmatS=optmat;
%             optmatS(a,:)=[];
%             AICminvecS=AICminvec;
%             AICminvecS(a)=[];
%             cdfRmatS=cdfRmat;
%             cdfRmatS(a,:)=[];
%             meanmatnormS=meanmatnorm;
%             meanmatnormS(a,:)=[];
%         else
%             %disp('fsdfsdf')
%         end
%     else
%         concvecS=concvec;
%         trymatoptrsgridS=trymatoptrsgrid;
%         optmatS=optmat; 
%         AICminvecS=AICminvec;
%         cdfRmatS=cdfRmat; 
%         meanmatnormS=meanmatnorm;
%     end
% end

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
    plot(trymatoptrsgridS(a,1:AICminvecS(a)),optmatS(a,1:AICminvecS(a)),'LineWidth',2,'Color',clr(a,:))
end
xlabel('Sensitivity to Treatment s')
ylabel('Recovered Proportion of Population')
legend(Legend, 'FontSize',12)
set(gca,"FontSize",20)


% %PLOT CDFS CORRECTLY
% figure
% for a=1:length(concvecS)
%     hold on
%     plot(trymatoptrsgridS(a,1:AICminvecS(a)),cdfRmatS(a,1:AICminvecS(a)),'LineWidth',2,'Color',clr(a,:))
% end
% xlabel('Sensitivity to Treatment s')
% ylabel('Cumulative Recovered Prop. of Pop.')
% legend(Legend, 'FontSize',12,'Location','northwest')
% set(gca,"FontSize",20)
% 
%plot cdfs again w stair fn
figure
for a=1:length(concvecS)
    hold on
    %X = linspace(0,4*pi,40);
    clr = cool(length(concvecS));
    X=trymatoptrsgridS(a,1:AICminvecS(a));
    %Y = sin(X);
    Y=cdfRmatS(a,1:AICminvecS(a));
    stairs(X,Y,'Color',clr(a,:),'LineWidth',2)
    ylim([0 1])
    xlabel('Sensitivity to Treatment s')
ylabel('Cumulative Recovered Prop. of Pop.')
legend(Legend, 'FontSize',12,'Location','northwest')
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
    data=meanmatnorm(a,1:tspanvec(a))
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
    clr = cool(length(concvecS));    %# LINES colormap
    plot(tspanS,data,'*','LineWidth',2,'Color',clr(a,:),'HandleVisibility', 'off')
    hold on
    plot(tspanS,weightedsol,'LineWidth',2,'Color',clr(a,:))
    hold on
    xlabel('Time (hours)')
    ylabel('Tumor Growth Data and Model Fit')
    legend(Legend, 'FontSize',12,'Location','northwest')
    set(gca,"FontSize",20)
end

%plot just first concentration data with forward solution to sus out issues
%with rhomax, and then do smth similar for rhomin. 
% if choosedata =='R500'
%  R500C1sol=[0.0320
%     0.0390
%     0.0474
%     0.0576
%     0.0697
%     0.0842
%     0.1014
%     0.1216
%     0.1452
%     0.1725
%     0.2037
%     0.2389
%     0.2780
%     0.3209];
% 
%  R500C1sol=[0.032011761358247
%    0.038995149857732
%    0.047427335853808
%    0.057573632097170
%    0.069731646596216
%    0.084227650322267
%    0.101408643141092
%    0.121628793035960
%    0.145229057846475
%    0.172509261297746
%    0.203692722723410
%    0.238885744993538
%    0.278036179211446
%    0.320897692527376];
%  figure
%  plot(tspan,R500C1sol)
%  hold on
%  plot(tspan,wsolmat(1,:))
%  hold on
%  plot(tspan,meanmatnorm(1,:),'*')
%  legend('from find max rho','from main code w optpar','data')
% 
%  maxrhoerr=R500C1sol' - wsolmat(1,:)
% 
% elseif choosedata == 'R250'
%  %    R250C11Sol=[0.026424453248301
%  %   0.028345775357960
%  %   0.030402434365977
%  %   0.032603309442530
%  %   0.034957764639875
%  %   0.037475660662994
%  %   0.040167364294590
%  %   0.043043754838706
%  %   0.046116227566958
%  %   0.049396693356351
%  %   0.052897573582812
%  %   0.056631789491027
%  %   0.060612745491865
%  %   0.064854305840611]
%  % 
%  %    figure
%  % plot(tspan,R250C11Sol)
%  % hold on
%  % plot(tspan,wsolmat(11,:))
%  % hold on
%  % plot(tspan,meanmatnorm(a,:),'*')
%  % legend('from finding min r','from foward prob w optpar','data')
%  % ylabel('tumor volume')
%  % xlabel('time')
%  % 
%  % minrhoerr=R250C11Sol' - wsolmat(11,:)
%  disp('bleh')
% elseif choosedata == 'S100'
% 
%   S100C10sol =[0.008203872244304
%    0.008088882040850
%    0.007975503606384
%    0.007863714349424
%    0.007753491995149
%    0.007644814580948
%    0.007537660452037
%    0.007432008257174
%    0.007327836944393
%    0.007225125756800
%    0.007123854228430
%    0.007024002180185
%    0.006925549715820
%    0.006828477217971];
%   figure
%  plot(tspan,S100C10sol)
%  hold on
%  plot(tspan,wsolmat(10,:))
%  hold on
%  plot(tspan,meanmatnorm(10,:))
%  legend('soln w min rho','sol from main code','data')
% 
% 
% else
%     disp('bleh')
% end

%NEW IDEA: SCALE K BY DIFFERENCE BETWEEN LARGEST AND SMALLEST POP OR SMTH?
%make kvec more linear by using log?



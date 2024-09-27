%Run inverse problem individually on all replicates from a given dataset 

clc
clear all

choosedata=['R500']; %choose between R250, R500, S500, and S1000
%or BF11,BF12,BF21,BF41
kmax=0.0077;
kmin=0;
rhomax=0.0682;
rhomin=0.0375;

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
else
    disp('??')
end
dosevec= [0, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];
dosevecshort= [0.00, 0.03,0.05, 0.12, 0.22, 0.34,0.46,1.19,2.44,3.71,5.01];
dosevecstring= ["0.00", "0.03","0.05", "0.12", "0.22", "0.34","0.46","1.19","2.44","3.71","5.01"];

kvec=kmax.*(dosevec./dosevec(11));

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

load('MONOCLONAL_DATA.mat');
load('Mixture_Data.mat');

if choosedata=='R250'
    T=RESISTANT_250_BF;
elseif choosedata=='R500'
    T=RESISTANT_500_BF;
elseif choosedata=='S500'
    T=SENSITIVE_500_BF;
elseif choosedata=='S100'
    T=SENSITIVE_1000_BF;
elseif choosedata=='BF11'
    T=BF_11;
elseif choosedata=='BF12'
    T=BF_12;
elseif choosedata=='BF21'
    T=BF_21;
elseif choosedata=='BF41'
    T=BF_41;
else
    disp('??')
end

sz = size(T);
tsize=sz(3);
concsize=sz(2);
tvec=linspace(1,tsize,tsize);
concvec=linspace(1,concsize,concsize);
%T is (replicate, concentration, time)


%choose dosage:
dosage = 11;


tinit=9;
tfinal=48;
tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

if choosedata=='R250'
    maxval=3*max(T,[],'all');
elseif choosedata=='R500'
    maxval=3*max(T,[],'all');
elseif choosedata=='S500'
    maxval=10*max(T,[],'all');
elseif choosedata=='S100'
    maxval=10*max(T,[],'all');
else 
    maxval=5*max(T,[],'all');
end

%set recovered mesh to four points 
rpoints = 4;

clr = hsv(sz(1));

% figure
% datavec = squeeze(T(4,dosage,:));
% tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);
% datavec = datavec/maxval;
% %plot replicates together
% array = isnan(datavec);
% datavecnew = datavec(~array);
% tspannew = tspan(~array);
% subplot(1,2,1)
% plot(tspannew, datavecnew,'o','Color',clr(4,:))
% hold on
% [gls_optpar,converge_flag,AIC_GLS,weightedsol,rsgrid,finalerr]=GLSInverseFnNData(rhovec(dosage),kvec(dosage),datavecnew(1),tspannew,rpoints,datavecnew');
% plot(tspannew,weightedsol,'Color',clr(4,:))
%     hold on
%     subplot(1,2,2)
%     plot(rsgrid, gls_optpar,'Color',clr(4,:))
%     hold on

    %% 
%for loop of each replicate
figure
for i=1:sz(1)
    datavec = squeeze(T(i,dosage,:));
    tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);
    %get rid of nans
    %normalize the data
    datavec = datavec/maxval;
    %plot replicates together
    array = isnan(datavec);
    datavecnew = datavec(~array);
    tspannew = tspan(~array);
    subplot(1,2,1)
    plot(tspannew, datavecnew,'o','Color',clr(i,:))
    hold on
    [gls_optpar,converge_flag,AIC_GLS,weightedsol,rsgrid,finalerr]=GLSInverseFnNData(rhovec(dosage),kvec(dosage),datavecnew(1),tspannew,rpoints,datavecnew');
    plot(tspannew,weightedsol,'Color',clr(i,:))
    hold on
    subplot(1,2,2)
    plot(rsgrid, gls_optpar,'Color',clr(i,:))
    hold on
end
%% 
meanvec = zeros(tsize,1);
meanmat = zeros(concsize,tsize);
concnum=length(concvec);

for i=tvec
    for j=concvec
        m=mean(T(:,j,i),'omitnan');
    %m = mean(T(i,:),'omitnan');
    %meanvec(i)=m;
        meanmat(j,i)=m;
    end
end


if choosedata=='R250'
    maxval=3*max(meanmat,[],'all');
elseif choosedata=='R500'
    maxval=3*max(meanmat,[],'all');
elseif choosedata=='S500'
    maxval=10*max(meanmat,[],'all');
elseif choosedata=='S100'
    maxval=10*max(meanmat,[],'all');
else 
    maxval=5*max(meanmat,[],'all');
end

%normalize the data
meanmatnorm=zeros(size(meanmat));
for a=concvec
    for t=tvec
        s=meanmat(a,t)/maxval; %multiply by 10 because data still growing?
        meanmatnorm(a,t)=s;
    end 
end

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
else
    disp('??')
end
dosevec= [0, 0.0260,0.0536, 0.1094, 0.2244, 0.3410,0.4603,1.1908,2.4430,3.7133,5.0119];
dosevecshort= [0.00, 0.03,0.05, 0.12, 0.22, 0.34,0.46,1.19,2.44,3.71,5.01];
dosevecstring= ["0.00", "0.03","0.05", "0.12", "0.22", "0.34","0.46","1.19","2.44","3.71","5.01"];

kvec=kmax.*(dosevec./dosevec(11));

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

%[meanmatnorm,concvec]=ScaleData(choosedata);
%need to scale data according to 3 or 10 or 5 times the max val


%[meanmatnorm,isnanmat, tspanvec, tspanmat,concvecA]=DropNaNsFn(concvec,meanmatnorm,tspan);
%need to drop nans from the data 
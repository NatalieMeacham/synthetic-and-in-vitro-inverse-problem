function [meanmatnorm,concvec]=ScaleData(choosedata)

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
%concvec=linspace(2,concsize,concsize-1);
%concvec=linspace(1,concsize-1,concsize-1)
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
    maxval=10*max(meanmat,[],'all'); %normally 10
elseif choosedata=='S100'
    maxval=10*max(meanmat,[],'all'); %normally 10
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



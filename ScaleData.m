function [meanmatnorm,concvec,errmat]=ScaleData(choosedata)

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
meanvec = zeros(tsize,1);
meanmat = zeros(concsize,tsize);
concnum=length(concvec);
sdmat = zeros(concsize,tsize);
errmat = zeros(concsize,tsize);

for i=tvec
    for j=concvec
        [s,m] = std(T(:,j,i),'omitnan');
        meanmat(j,i)= m;
        sdmat(j,i) = s;
        errmat(j,i) = s/sqrt(length(T(:,j,i)) - sum(isnan(T(:,j,i))));
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

maxval1=max(RESISTANT_250_BF);
maxval2=max(RESISTANT_500_BF);
maxval3=max(SENSITIVE_500_BF);
maxval4=max(SENSITIVE_1000_BF);
maxlist=[maxval1,maxval2, maxval3,maxval4];
maxval=10*max(maxlist,[],'all');

meanmatnorm=zeros(size(meanmat));

Tnormalized = T./maxval;

for a=concvec
    for t=tvec
        [s,m] = std(Tnormalized(:,a,t),'omitnan');
        meanmatnorm(a,t)= m;
        sdmat(a,t) = s;
        errmat(a,t) = s/sqrt(length(T(:,a,t)) - sum(isnan(T(:,a,t))));
    end 
end



clc
clear all
points=150; %this used to be 50
disttype = 'Bigaussian';
rho=1; %make these inputs
k=1.5; %maximal death rate due to treatment %bigger than rho
y0=.2; %note c is normalized so needs to start between 0 and 1
tfinal=10;
tpoints=25;
noisesize=0.0; %keep it w 0 noise for now
tspan=linspace(0,tfinal,tpoints);
rpoints = 15; %ER recommended like 15
a=0;
b=1;
sgrid=linspace(0,1,points);

[sprobs] = DistFn(disttype,sgrid,a,b)

[t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan)

propdata=weightedsol.*(1 + noisesize*(-1 + (2).*rand(size(weightedsol))));

[gls_optpar,converge_flag,AIC_GLS,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints,propdata)

%next step is to make pdf, new cdf, and fit plots for a few different noise
%sizes 
%also need to remake other fig (the one w the different dists)
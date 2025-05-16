clc
clear all

load('MONOCLONAL_DATA.mat');
choosedata='S500';

tinit=9;
tfinal=48;
tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

%choose concentration
chooseconc=10; 

[meanmatnorm,concvec]=ScaleData(choosedata);

[meanmatnorm,isnanmat, tspanvec, tspanmat]=DropNaNsFn(concvec,meanmatnorm,tspan);

data=meanmatnorm(chooseconc,1:tspanvec(chooseconc));
tspan=tspanmat(chooseconc,1:tspanvec(chooseconc));
data;
tspan;

IC=data(1);
Aeq=[]; 
beq=[];
lb=0; 
ub=1;

k0=0.1;

gammas=1;
weights=ones(size(tspan));

gls_optpar=k0; %initial optpar 
old_gls_optpar=k0; %initial old optpar 
tol=1e-4; %tolerance: keeps weights above a certain size  
maxits=2000;  
minits=10; 
partol=0.1; %another tolerance, convergence criteria  
parchange=100; %initialize parameter change, start w smth large  
oldparchange=100; %initialize old parameter change  
ii=1; %initialize number of iterations 

while ii<maxits & parchange > partol & oldparchange > partol | ii< minits 
errtomin=@(par)ErrorFnFindK(par, data,tspan,IC,weights);
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); 
[gls_optpar,~,~,~]=fmincon(errtomin, k0, [], [], Aeq, beq, lb, ub, [], options); 
converge_flag=1;
[~,weights]=ForwardFnFindK(gls_optpar,tspan,IC);
weights=weights';
weights(weights<tol)=0; 
weights(weights>tol)=weights(weights>tol).^(-2*gammas); 
inds=old_gls_optpar>1e-10; 
weights=full(weights); 
parchange =1/(2)*sum((abs(gls_optpar(inds)-old_gls_optpar(inds))./old_gls_optpar(inds))); 
ii = ii+1; 
old_gls_optpar=gls_optpar; 
end 
disp('GLS Estimation') 
gls_optpar

figure
[t,c]=ForwardFnFindK(gls_optpar,tspan,IC);
plot(tspan,data,'*','LineWidth',2)
hold on
plot(tspan,c,'LineWidth',2)
legend('Data','Soln. w/ Max. k')
xlabel('Time (Hours)')
ylabel('Most Sensitive S. Pop.')
set(gca,"FontSize",20)

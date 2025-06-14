function[gls_optpar,converge_flag,AIC_GLS,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints,propdata)

%Function to run the inverse problem for a given dataset where the mesh of
%the recovered senstivity distribution is given

%INPUTS
    %points: number of points in original sensitivity distribution
    %disttype: type of initial distribution
    %rho: growth rate
    %k: death rate
    %y0: initial tumor volume
    %tfinal: final time point
    %tpoints: number of time points
    %noisesize: size of noise applied to tumor growth curve to get data
    %rpoints: number of points in recovered distribution
    %propdata: proportional error data (synthetic data)

%OUTPUTS
    %gls_optpar: recovered sensitivity distribution
    %converge_flag: convergence flag 
    %AIC_GLS: AIC score for given gls_optpar
    %sgrid: mesh for original sensitivity distribution
    %sprobs: original sensitivity distribution
    %weightedsol: aggregated tumor volume from original distribution
    %rsgrid: mesh for recovered distribution 

%FORWARD PROBLEM
%Create synthetic data
 a=0; 
 b=1;
 sgrid=linspace(a,b,points);
 
%Create original distribution
[sprobs] = DistFn2(disttype,sgrid,a,b);

%Compute tspan for forward solver to use
tspan=linspace(0,tfinal,tpoints); 

%Run forward function 
[t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);

%Set up mesh for recovered distribution
rsgrid = linspace(a,b,rpoints);

%INVERSE PROBLEM
%Set up constraints for constrained optimization (recovered s dist. must sum to 1 and have nonneg vals)
%Must sum to 1:
Aeq=zeros(1,rpoints); 
Aeq(1,:)=1;
beq=ones(1,1);
%Only nonneg vals:
lb=zeros(rpoints,1); 
ub=zeros(rpoints,1);
ub(1:rpoints)=1;

%Make normalized uniform dist. initial guess for s0
s0=ones(rpoints,1);
s0=s0/sum(s0); 
init_guess=s0';

%Weights for use in PEerrorfn 
gammas=1;
weights=ones(size(weightedsol)); 

gls_optpar=init_guess; %initial optpar 
old_gls_optpar=init_guess; %initial old optpar 
tol=1e-4; %tolerance: keeps weights above a certain size  
maxits=2000;  
minits=10; 
partol=0.1; %another tolerance, convergence criteria  
parchange=100; %initialize parameter change, start w smth large  
oldparchange=100; %initialize old parameter change  
ii=1; %initialize number of iterations  

%While loop within tolerances to find optimal recovered distribution
while ii<maxits & parchange > partol & oldparchange > partol | ii< minits 
[t, gencmatC,~] = ForwardFunctionN(rsgrid, ones(size(rsgrid)), rho, k, y0, tspan);
errtomin=@(par)PEerrorfn(par,gencmatC,propdata,weights);
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); 
[gls_optpar,~,converge_flag,~]=fmincon(errtomin, s0, [], [], Aeq, beq, lb, ub, [], options);
[~, ~,weights] = ForwardFunctionN(rsgrid, gls_optpar', rho, k, y0, tspan); 
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

%AIC
%final error between data and aggregated cell density from recovered mesh
[finalerr,~]=PEerrorfn(gls_optpar, gencmatC, propdata,weights);

%Compute AIC score for the given points and rpoints here:
AIC_GLS = (tpoints)*log(sum((weights).*((propdata-(gencmatC*gls_optpar)').^2))/tpoints) + 2*(rpoints + 1);

end
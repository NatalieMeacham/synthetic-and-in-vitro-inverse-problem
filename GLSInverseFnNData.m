function[gls_optpar,converge_flag,AIC_GLS,weightedsol,rsgrid,finalerr]=GLSInverseFnNData(rho,k,y0,tspan,rpoints,data)

% Define sensitivity grid
a=0; 
b=1;
rsgrid=linspace(a,b,rpoints);


%Set up constraints for constrained optimization (ie s must sum to 1 and have nonneg vals)
%Note that these are constraints for what we're trying to recover
%Must sum to 1:
Aeq=zeros(1,rpoints); %this chunk does not want to be fn inputs 
Aeq(1,:)=1;
beq=ones(1,1);
%Only nonneg vals:
lb=zeros(rpoints,1); 
ub=zeros(rpoints,1);
ub(1:rpoints)=1;

%Define initial s which is uniform s0
s0=ones(rpoints,1); %to get a uniform dist.
s0=s0/sum(s0); %normalize
init_guess=s0';

%Weights for use in PEerrorfn 
gammas=1;
weights=ones(size(tspan)); %doesn't need to be input

gls_optpar=init_guess; %initial optpar 
old_gls_optpar=init_guess; %initial old optpar 
tol=1e-4; %tolerance: keeps weights above a certain size  
maxits=2000;  
minits=10; 
partol=0.1; %another tolerance, convergence criteria  
parchange=100; %initialize parameter change, start w smth large  
oldparchange=100; %initialize old parameter change  
ii=1; %initialize number of iterations  

while ii<maxits & parchange > partol & oldparchange > partol | ii< minits 

%gls_error_estimate = @(par)gls_formulation(par,prop_data,ts,y0,weights); 
[t, gencmatC,weightedsol] = ForwardFunctionN(rsgrid, ones(size(rsgrid)), rho, k, y0, tspan); %WHY ONES? because we're only using the genmatC and don't care about the the weights? 
errtomin=@(par)PEerrorfn(par,gencmatC,data,weights);
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); 
[gls_optpar,~,~,~]=fmincon(errtomin, s0, [], [], Aeq, beq, lb, ub, [], options); %without converge flag
converge_flag=1;
%pause
[~, ~,weights] = ForwardFunctionN(rsgrid, gls_optpar', rho, k, y0, tspan); 
%weights is from the weightedsol so here would be c
weights(weights<tol)=0; 
weights(weights>tol)=weights(weights>tol).^(-2*gammas); %note that the weights already get taken to the -2 here
inds=old_gls_optpar>1e-10; 
weights=full(weights); 
parchange =1/(2)*sum((abs(gls_optpar(inds)-old_gls_optpar(inds))./old_gls_optpar(inds))); 
ii = ii+1; 
old_gls_optpar=gls_optpar; 
oldparchange = parchange;
end 
disp('GLS Estimation') 
gls_optpar; 
 
%AIC SECTION
%take error where s is given by optweight for AIC
%[current_err,~] = errorfunc_discrete_sparse(optweight,fullsol,agg_sol,weights);
[finalerr,~]=PEerrorfn(gls_optpar, gencmatC, data,weights);

tpoints=length(tspan); 

%Compute AIC score for the given points and rpoints here:
AIC_GLS = (tpoints)*log(finalerr/tpoints) + 2*(rpoints + 1);

end
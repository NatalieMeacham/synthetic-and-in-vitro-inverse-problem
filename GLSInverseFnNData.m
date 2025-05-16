function[gls_optpar,converge_flag,AIC_GLS,weightedsol,rsgrid,finalerr]=GLSInverseFnNData(rho,k,y0,tspan,rpoints,data)

a=0; 
b=1;
rsgrid=linspace(a,b,rpoints);

%must sum to 1:
Aeq=zeros(1,rpoints); %this chunk does not want to be fn inputs 
Aeq(1,:)=1;
beq=ones(1,1);
%only nonneg vals:
lb=zeros(rpoints,1); 
ub=zeros(rpoints,1);
ub(1:rpoints)=1;

%uniform initial guess s0
s0=ones(rpoints,1); %to get a uniform dist.
s0=s0/sum(s0); %normalize
init_guess=s0;

%weights for use in PEerrorfn 
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
[t, gencmatC,weightedsol] = ForwardFunctionN(rsgrid, ones(size(rsgrid)), rho, k, y0, tspan);  
errtomin=@(par)PEerrorfn(par,gencmatC,data,weights);
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); 
[gls_optpar,~,~,~]=fmincon(errtomin, old_gls_optpar, [], [], Aeq, beq, lb, ub, [], options); 
converge_flag=1;
[~, ~,weights] = ForwardFunctionN(rsgrid, gls_optpar', rho, k, y0, tspan); 
weights(weights<tol)=0; 
weights(weights>tol)=weights(weights>tol).^(-2*gammas); 
inds=old_gls_optpar>1e-10; 
weights=full(weights); 
parchange =1/(2)*sum((abs(gls_optpar(inds)-old_gls_optpar(inds))./old_gls_optpar(inds))); 
ii = ii+1; 
old_gls_optpar=gls_optpar; 
oldparchange = parchange;
end 

disp('GLS Estimation') 
gls_optpar; 
 
%compute error from final estimate and then compute AIC
[finalerr,weightedsol]=PEerrorfn(gls_optpar, gencmatC, data,weights);
tpoints=length(tspan); 
AIC_GLS = (tpoints)*log(finalerr/tpoints) + 2*(rpoints + 1);

end
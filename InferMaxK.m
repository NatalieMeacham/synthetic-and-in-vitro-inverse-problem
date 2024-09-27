clc
clear all

load('MONOCLONAL_DATA.mat');
%T=RESISTANT_250_BF; %first row has NaNs, use concvec=linspace(2,concsize,concsize-1);
%T=RESISTANT_500_BF; %last row has NaN, use concvec=linspace(1,concsize-1,concsize-1);
%T=SENSITIVE_500_BF; %no NaN
%T=SENSITIVE_1000_BF; %last row has NaN, use concvec=linspace(1,concsize-1,concsize-1);
choosedata='S500';
% constdata=T;
% sz = size(T);
% tsize=sz(3);
% concsize=sz(2);
% tvec=linspace(1,tsize,tsize);
% concvec=linspace(1,concsize,concsize);
% meanvec = zeros(tsize,1);
% meanmat = ones(concsize,tsize);
tinit=9;
tfinal=48;
tspan=linspace(tinit,tfinal,(tfinal-tinit)/3 + 1);

% concnum=length(concvec);

chooseconc=10;

% for i=tvec
%     for j=concvec
%         m=mean(T(:,j,i),'omitnan');
%     %m = mean(T(i,:),'omitnan');
%     %meanvec(i)=m;
%         meanmat(j,i)=m;
%     end
% end

[meanmatnorm,concvec]=ScaleData(choosedata);
meanmatnorm
%keyboard
%meanmatnorm
[meanmatnorm,isnanmat, tspanvec, tspanmat]=DropNaNsFn(concvec,meanmatnorm,tspan);
%keyboard
meanmatnorm

% %HERE WE ARE TRYING SMTH DIFFERENT W A INSTEAD OF CONCNUM
% maxval = 10*max(meanmat,[],'all');
% meanmatnorm=zeros(size(meanmat));
% for a=concvec
%     for t=tvec
%         s=meanmat(a,t)/maxval; %multiply by 10 because data still growing?
%         meanmatnorm(a,t)=s;
%     end
% 
% end
% 
%replace any NaNs with 1s
% TFvec=zeros(length(concvec),1);
% for a=concvec
%     TFvec(a)=anynan(meanmatnorm(a,:))
%     %TFvec=[TFvec, TF]
% end
% % 
% for a=concvec
%     if TFvec(a)==1
%         meanmatnorm(a,:)=1;
%         %disp('updated meanmatnorm with ones')
%     else 
%         %disp('chill')
%     end
% end

% [meanmatnorm,concvec]=ScaleData(choosedata);
% meanmatnorm
% keyboard
% %meanmatnorm
% [meanmatnorm,isnanmat, tspanvec, tspanmat]=DropNaNsFn(concvec,meanmatnorm,tspan)
% meanmatnorm
% data=meanmatnorm(chooseconc,:);
% data
% keyboard
data=meanmatnorm(chooseconc,1:tspanvec(chooseconc));
tspan=tspanmat(chooseconc,1:tspanvec(chooseconc));
data
tspan
%keyboard
%keyboard

IC=data(1);
Aeq=[]; %this chunk does not want to be fn inputs 
%Aeq(1,:)=[];
beq=[];
%Only nonneg vals:
lb=0; 
ub=1;
%ub(1:rpoints)=1;

k0=0.0005;

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
%gls_error_estimate = @(par)gls_formulation(par,prop_data,ts,y0,weights); 
%[t, gencmatC,weightedsol] = ForwardFunctionN(rsgrid, ones(size(rsgrid)), rho, k, y0, tspan);
%[t,c]=ForwardFnFindK(1,tspan,IC);
%errtomin=@(par)PEerrorfn(par,gencmatC,data,weights);
%errtomin=@(par)ErrFn2K(par, data,tspan,IC,weights);
errtomin=@(par)ErrFn2K(par, data,tspan,IC,weights);
options = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIterations',550); 
%gls_optpar=FMINCON (set up Aeq etc.) 
%[gls_optpar,~,converge_flag,~]=fmincon(errtomin, s0, [], [], Aeq, beq, lb,
%ub, [], options); %with converge flag
weights;
%keyboard
[gls_optpar,~,~,~]=fmincon(errtomin, k0, [], [], Aeq, beq, lb, ub, [], options); %without converge flag
gls_optpar
%keyboard
converge_flag=1;
%pause
%[~, ~,weights] = ForwardFunctionN(rsgrid, gls_optpar', rho, k, y0, tspan);
[~,weights]=ForwardFnFindK(gls_optpar,tspan,IC)
weights=weights';
%keyboard
%weights is from the weightedsol so here would be c
weights(weights<tol)=0; 
weights(weights>tol)=weights(weights>tol).^(-2*gammas); %note that the weights already get taken to the -2 here
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
xlabel('Time')
ylabel('Most Sensitive S. Population')
set(gca,"FontSize",20)

clear all
points=100;
disttype='TwoPoints';
rho=0.3;
k=0.45;
y0=0.2;
tfinal=10;
tpoints=25;
noisesize=0;
rpoints=10;
a=0;
b=1;
tspan=linspace(0,tfinal,tpoints);
sgrid=linspace(0,1,points);
sprobs = DistFn2(disttype,sgrid,a,b);


%create propdata
[t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan);
[t, cmatR,~] = RK4FunctionC(linspace(0,1,rpoints), ones(1,rpoints), rho, k, y0, tspan);
propdata=weightedsol.*(1 + noisesize*randn(size(weightedsol)));

[gls_optpar,converge_flag,AIC_GLS,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnC(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints,propdata)

figure
yyaxis left
stem(rsgrid,gls_optpar,'--','MarkerSize',8,'LineWidth',2,'Color','blue') %,'MarkerSize,'8,'LineWidth',2,'Color,''blue'))
ylabel('Recovered Prop. of Pop.')
hold on
yyaxis right
plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
xlabel('Sensitivity to Treatment')
ylabel('Original Prop. of Pop.')
set(gca,"FontSize",20)

%cdf
if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') == 1 
    cdfR=cumsum(gls_optpar)
    cdfS=cumsum(sprobs)
else
    cdfR=cumsum(gls_optpar)
    cdfS=cumtrapz(sgrid,sprobs)
end
figure
stairs(rsgrid,cdfR,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
hold on
stairs(sgrid,cdfS,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
legend('original CDF','recovered CDF')
set(gca,"FontSize",20)
xlabel('Sensitivity to Treatment')
ylabel('Cumulative Prop. of Pop.')

figure
plot(t,propdata,'r*','MarkerSize',10,'LineWidth',1.5)
hold on
plot(t,cmatR*gls_optpar,'b','LineWidth',2)
legend('propdata','recovered fit')
xlabel('Time')
ylabel('Aggregated Tumor Volume')
set(gca,"FontSize",20)
%red original, blue recovered 
clc
clear all
points=101; %this used to be 50
disttype = ['TwoPoints'];
% rho=1; %make these inputs 
% k=1.5; %maximal death rate due to treatment %bigger than rho
% y0=.2; %note c is normalized so needs to start between 0 and 1
tfinal=50;
tpoints=100; %typically 100 for paper
noisesize=0.05; %keep it w 0 noise for now
tspan=linspace(0,tfinal,tpoints);
rpoints = 12; 
pointsstr=string(points);
noisestr=string(noisesize);

rho=0.3;
k=0.45;
y0=.2; 
% a=0;
% b=1;
% sgrid=linspace(a,b,points);
% [sprobs] = DistFn(disttype,sgrid,a,b);
% 
% [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);
% 
% rsgrid=linspace(0,1,rpoints);
% s0=ones(length(rsgrid),1);
% s0=s0/sum(s0);
% s0=s0';
% [t, cmat,weightedsol] = ForwardFunctionN(rsgrid, s0, rho, k, y0, tspan);

a=0; 
b=1;
sgrid=linspace(a,b,points);

%Create original distribution
[sprobs] = DistFn2(disttype,sgrid,a,b);
% sprobs
% keyboard

[t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);
%propdata= weightedsol.*(1+noisesize*rand(size(weightedsol))); 
%propdata=weightedsol.*(1 + noisesize*(-1 + (2).*rand(size(weightedsol))));
propdata=weightedsol.*(1 + noisesize*randn(size(weightedsol)));

% figure
% plot(t,propdata)

% figure
% plot(tspan,weightedsol)
% hold on
% plot(tspan,propdata,'*')
% keyboard
%%
%[gls_optpar,converge_flag,AIC_GLS,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints)
[gls_optpar,converge_flag,AIC_GLS,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints,propdata)

% figure
% plot(t,weightedsol)
%%

%copied from OLSInverseScript
rpointsvec=4:1:30;
AICvec=zeros(length(rpointsvec),1);
%ii=100:50:1000
for i=1:length(rpointsvec)
    %[optweight,~,AIC_OLS,~,~,~,rsgrid]=OLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpointsvec(i));
    %[optweight,~,AIC_GLS,~,~,~,rsgrid]=GLSInverseFn(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpointsvec(i));
    [gls_optpar,~,AIC_GLS,~,~,~,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpointsvec(i),propdata);
    AICvec(i)=AIC_GLS;
end

[M,I]=min(AICvec);
optrpoints=rpointsvec(I);

%get optweight from OLSInverseFn for optrpoints
%SHOULD THIS BE OPTRPOINTS NOT RPOINTS???
%[optweightfromAIC,~,~,sgrid,sprobs,weightedsol,rsgrid]=OLSInverseFn(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpoints);
%[optweightfromAIC,~,~,sgrid,sprobs,weightedsol,rsgrid]=OLSInverseFn(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,optrpoints);
%[optweightfromAIC,~,~,sgrid,sprobs,weightedsol,rsgrid]=OLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,optrpoints);
%[optweightfromAIC,~,~,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFn(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,optrpoints)
[optweightfromAIC,~,~,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,optrpoints,propdata);

figure
plot(rpointsvec,AICvec,'LineWidth',2,'Color','blue')
xlabel('Number of Points in Recovered Dist.')
ylabel('AIC Score')
%title('AIC Score for Differing Meshes')
set(gca,"FontSize",20)
hold on
plot(rpointsvec(I),M,'m*','LineWidth',2,'MarkerSize',20)
hold off
legend('AIC','Min. AIC','Location','northeast')
AICfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'AIC','.jpg');
saveas(gcf,AICfiglabel);


%comparison of the og and recovered dists 

% %compare output of fmincon with original dist
% figure
% yyaxis left
% %plot(rsgrid,optweightfromAIC,'-*','LineWidth',2)
% plot(rsgrid,optweightfromAIC,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
% %ylabel('Recovered Prop. of Pop.') %original
% ylabel('Recovered Proportion of Population')
% ylimval=max(max(sprobs),max(optweightfromAIC));
% if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1 
%     ylim([0 ylimval])
% end
% if strcmp(disttype,'Uniform') == 1
%     ylim([0 4*median(optweightfromAIC)])
% end
% hold on
% yyaxis right
% plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
% legend('Recovered','Original','Location','northeast')
% xlabel('Sensitivity to Treatment {\it s}')
% %ylabel('Original Prop. of Pop.') %original version
% ylabel('Original Proportion of Population')
% if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
%     ylim([0 ylimval])
% end
% if strcmp(disttype,'Uniform') == 1
%     ylim([0 4*median(sprobs)])
% end
% set(gca,"FontSize",20)
% PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists','.jpg');
% saveas(gcf,PMFfiglabel);
%title('Recovered vs. Original Distributions')
%figurename = strcat('F','rp',num2str(rpoints),'sp',num2str(points));
%saveas(gcf,figurename)


%compare output of fmincon with original dist (no lines between recovered pts)
figure
yyaxis left
%plot(rsgrid,optweightfromAIC,'-*','LineWidth',2)
%plot(rsgrid,optweightfromAIC,'o','MarkerSize',8,'LineWidth',2,'Color','blue')
stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) %,'MarkerSize,'8,'LineWidth',2,'Color,''blue')
%ylabel('Recovered Prop. of Pop.') %original
ylabel('Recovered Proportion of Population')
ylimval=max(max(sprobs),max(optweightfromAIC));
if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1 
    ylim([0 ylimval])
end
if strcmp(disttype,'Uniform') == 1
    ylim([0 2.5*median(optweightfromAIC)])
end
hold on
yyaxis right
plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
legend('Recovered','Original','Location','northeast')
xlabel('Sensitivity to Treatment {\it s}')
%ylabel('Original Prop. of Pop.') %original version
ylabel('Original Proportion of Population')
if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
    ylim([0 ylimval])
end
if strcmp(disttype,'Uniform') == 1
    ylim([0 2.5*median(sprobs)])
end
set(gca,"FontSize",20)
PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists','.jpg');
saveas(gcf,PMFfiglabel);
PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists','.fig');
saveas(gcf,PMFfiglabel);

%make and plot cdfs to accommodate cts dists:
%if disttype == ['Normal'] || disttype == ['Uniform'] || disttype == ['Bigaussian']
cdfR=zeros(optrpoints,1);
cdfS=zeros(points,1);
if strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Uniform') ==1 || strcmp(disttype,'Bigaussian') ==1 
    %stuff
    %cdfS=zeros(points,1);
    partialsumS=0;
    cdfS(1)=(1/points)*sprobs(1);
    for i=2:points
        %partialsumS=sum(sprobs(1:i));
        partialsumS=trapz(sgrid(1:i),sprobs(1:i));
        cdfS(i)=partialsumS;
    end
    %cdfR=zeros(optrpoints,1);
    partialsumR=0;
    for i=1:optrpoints
        partialsumR=sum(optweightfromAIC(1:i));
        cdfR(i)=partialsumR;
    end
    disp('Here is a cdf for a continuous dist.')
elseif strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
    %cdfS=zeros(points,1);
    partialsumS=0;
    for i=1:points
        partialsumS=sum(sprobs(1:i));
        cdfS(i)=partialsumS;
    end
    %cdfR=zeros(optrpoints,1);
    partialsumR=0;
    for i=1:optrpoints
        partialsumR=sum(optweightfromAIC(1:i));
        cdfR(i)=partialsumR;
    end
    disp('Here is a cdf for a discrete dist.')
end


% %make and plot cdfs
% cdfS=zeros(points,1);
% partialsumS=0;
% for i=1:points
%     partialsumS=sum(sprobs(1:i));
%     cdfS(i)=partialsumS;
% end
% 
% cdfR=zeros(optrpoints,1);
% partialsumR=0;
% for i=1:optrpoints
%     partialsumR=sum(optweightfromAIC(1:i));
%     cdfR(i)=partialsumR;
% end

figure
    %X = linspace(0,4*pi,40);
    %clr = cool(length(concvecS));
    X=rsgrid;
    %Y = sin(X);
    Y=cdfR;
    stairs(X,Y,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
    ylim([0 1])
    hold on
    plot(sgrid, cdfS,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
    xlabel('Sensitivity to Treatment {\it s}')
%ylabel('Cumulative Recovered Prop. of Pop.') %original version
ylabel('Cumulative Recovered Proportion')
%legend(Legend, 'FontSize',12,'Location','northwest')
legend('Recovered','Original','Location','Southeast')
set(gca,"FontSize",20)
CDFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'CDF','.jpg');
saveas(gcf,CDFfiglabel);
CDFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'CDF','.fig');
saveas(gcf,CDFfiglabel);

% figure
%     %X = linspace(0,4*pi,40);
%     %clr = cool(length(concvecS));
%     X=rsgrid;
%     Y=cdfR;
%     X1=sgrid;
%     %Y = sin(X);
%     Y1=cdfS;
%     stairs(X,Y,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
%     ylim([0 1])
%     hold on
%     stairs(X1,Y1,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
%     %plot(sgrid, cdfS,'LineWidth',2)
%     xlabel('Sensitivity to Treatment s')
% ylabel('Cumulative Recovered Prop. of Pop.')
% %legend(Legend, 'FontSize',12,'Location','northwest')
% legend('Recovered','Original','Location','Southeast')
% set(gca,"FontSize",20)

% %CDF figure adjusted for cts versions (and all versions)
% figure
% yyaxis left
% X=rsgrid;
%     Y=cdfR;
%     X1=sgrid;
%     %Y = sin(X);
%     Y1=cdfS;
%     stairs(X,Y,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
% %plot(rsgrid,optweightfromAIC,'-*','LineWidth',2)
% %plot(rsgrid,optweightfromAIC,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
% ylabel('Recovered CDF')
% hold on
% yyaxis right
% stairs(X1,Y1,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
%     %plot(sgrid, cdfS,'LineWidth',2)
%     xlabel('Sensitivity to Treatment {\it s}')
% ylabel('Original CDF')
% %legend(Legend, 'FontSize',12,'Location','northwest')
% legend('Recovered','Original','Location','Southeast')
% set(gca,"FontSize",20)
% %title('Recovered vs. Original Distributions')
% %figurename = strcat('F','rp',num2str(rpoints),'sp',num2str(points));
% %saveas(gcf,figurename)



% figure
% plot(rsgrid, cdfR,'LineWidth',2)
% hold on
% plot(sgrid, cdfS,'LineWidth',2)
% legend('Recovered','Original','Location','Southeast')
% xlabel('Sensitivity to Treatment')
% ylabel('Cumulative Proportion of Population')
% set(gca,"FontSize",20)
% hold off
% 

%constdata=weightedsol + noisesize*(-1 + (2).*rand(size(weightedsol))); 
%compare weighted csol from optweight with original weightedsol
[t,~,sweightedsol] = ForwardFunctionN(rsgrid, optweightfromAIC', rho, k, y0, tspan);
figure 
plot(t,sweightedsol,'LineWidth',2,'Color','blue')
hold on
%plot(t,weightedsol','rd','LineWidth', 2,'MarkerSize',8)
plot(t,propdata','rd','LineWidth', 2,'MarkerSize',8)
hold off
legend('Est. Tumor Volume','Simulated Data','Location','southeast')
xlabel('Time')
ylabel('Aggregated Tumor Volume')
%title('Weighted Sols for Original vs. Recovered Dists.')
set(gca,"FontSize",20)
Fitfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Fit','.jpg');
saveas(gcf,Fitfiglabel);
Fitfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Fit','.fig');
saveas(gcf,Fitfiglabel);


%compare output of fmincon with original dist on same scaled axes (no lines between recovered pts)
if strcmp(disttype,'Uniform') == 1 || strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Bigaussian') == 1
    figure
    stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) %,'MarkerSize,'8,'LineWidth',2,'Color,''blue')
    ylabel('Proportion of Population')
    hold on
    plot(sgrid,sprobs*length(sprobs)/(sum(sprobs)*length(optweightfromAIC)),'--o','MarkerSize',8,'LineWidth',2,'Color','red')
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[0 limsy(2)]);
    legend('Recovered','Original','Location','northeast')
    xlabel('Sensitivity to Treatment {\it s}')
    set(gca,"FontSize",20)
    PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists_samescale','.jpg');
    saveas(gcf,PMFfiglabel);
    PMFfiglabel=strcat(disttype,'P',pointsstr,'N',noisestr,'Dists_samescale','.fig');
    saveas(gcf,PMFfiglabel);
end

% figure
% %yyaxis left
% %plot(rsgrid,optweightfromAIC,'-*','LineWidth',2)
% %plot(rsgrid,optweightfromAIC,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
% %ylabel('Proportion of Pop.')
% %hold on
% %yyaxis right
% plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
% %legend('Recovered Dist.','Original Dist.','Location','northeast')
% xlabel('Sensitivity to Treatment (s)')
% ylabel('Original Prop. of Pop.')
% ylim([0 1])
% set(gca,"FontSize",20)

clc
clear all

%tic
a=0;
b=1;
points=11;
sgrid=linspace(a,b,points);

%[sprobs] = DistFn(['OnePoint'],sgrid,a,b);
[sprobs] = DistFn('Normal',sgrid,a,b);

%Other parameters for forward solution 
rho=1;
k=1.5;
y0=.2; %note c is normalized so needs to start between 0 and 1
tfinal=10;
tpoints=25;
tspan=linspace(0,tfinal,tpoints); %this is so we can initialize matrix with known number of rows 

%[t, cmat] = ForwardFunction(sgrid, sprobs, rho, k, y0, tspan)
[t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);

%toc
%DUMPED CODE
%Gaussian curve over sgrid: 
% sigma=1;
% mu=1;
% %sprobs=1/(sigma .* sqrt(2*pi)) .*exp(-(sgrid - mu)^2/(2.*sigma^2))
% sprobs=zeros(1,length(sgrid));
% for i=1:length(sgrid)
%    sprobs(i)=1/(sigma .* sqrt(2*pi)) .*exp(-(sgrid(i) - mu)^2/(2.*sigma^2));
% end
% sprobs

%plot s distribution
    figure
    plot(sgrid, sprobs,'*','LineWidth',2)
    %p.LineWidth = 2;
    set(gca,"FontSize",20)
    %title('Proportions of Different s Values')
    xlabel('Sensitivity to Treatment (s)')
    ylabel('Proportion of Population')
    ylim([0,1])
    %ylim([0,0.15])

    figure
    for i=1:length(sgrid)
        hold on
        plot(sgrid(i), sprobs(i),'*','LineWidth',2)
    end
    set(gca,"FontSize",20)
    %title('Proportions of Different s Values')
    xlabel('Sensitivity to Treatment (s)')
    ylabel('Proportion of Population')
    ylim([0,1])

    figure
    for i=1:length(sgrid)
        hold on
        plot(tspan, cmat(:,i),'LineWidth',2)
    end
    set(gca,"FontSize",20)
    %title('Proportions of Different s Values')
    xlabel('Sensitivity to Treatment (s)')
    ylabel('Proportion of Population')
    ylim([0,1])

    
    figure
    plot(tspan,weightedsol,'LineWidth',2)
    xlabel('Time')
    ylabel('Aggregated Tumor Volume')
    ylim([0,1])
    set(gca,"FontSize",20)

    figure
    scale=0.2;
    %data=weightedsol + rand(size(tspan))*scale*weightedsol;
    data=weightedsol.*(1+scale*rand(size(weightedsol)));
    plot(tspan,data,'*','LineWidth',2)
    xlabel('Time')
    ylabel('Aggregated Tumor Volume with Noise')
    ylim([0,1])
    set(gca,"FontSize",20)


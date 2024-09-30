clc
clear all

%tic
a=0;
b=1;
points=11;
sgrid=linspace(a,b,points);

%[sprobs] = DistFn(['OnePoint'],sgrid,a,b);
[sprobs] = DistFn2('Uniform',sgrid,a,b);

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
    % figure
    % plot(sgrid, sprobs,'*','LineWidth',2)
    % %p.LineWidth = 2;
    % set(gca,"FontSize",20)
    % %title('Proportions of Different s Values')
    % xlabel('Sensitivity to Treatment (s)')
    % ylabel('Proportion of Population')
    % ylim([0,1])
    % %ylim([0,0.15])

    figure
    for i=1:length(sgrid)
        hold on
        plot(sgrid(i), sprobs(i),'o','LineWidth',3,'MarkerSize',12)
    end
    set(gca,"FontSize",20)
    %title('Proportions of Different s Values')
    xlabel('Sensitivity to Treatment (s)')
    ylabel('Proportion of Population')
    %ylim([0,1])

    %plot all curves
    figure
    for i=1:length(sgrid)
        hold on
        plot(t, cmat(:,i),'-','LineWidth',2)
        Legend{i}=strcat('s=',' ', num2str(sgrid(i)));
    end
    set(gca,"FontSize",20)
    %title('Proportions of Different s Values')
    xlabel('Time')
    ylabel('Subpopulation Growth')
    legend(Legend)
    ylim([0,1])

    %plot weightedsol
    figure
    plot(t,weightedsol,'k','LineWidth',2)
    xlabel('Time')
    ylabel('Aggregated Volume')
    ylim([0,1])
    set(gca,"FontSize",20)

    % %plot a single solution like s=0.3
    % figure
    % plot(t,cmat(:,3),'r','LineWidth',2)
    % xlabel('Time')
    % ylabel('Subpopulation Growth')
    % legend('s=0.3')
    % ylim([0,1])
    % set(gca,"FontSize",20)

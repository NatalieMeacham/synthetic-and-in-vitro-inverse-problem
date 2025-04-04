a=0;
b=1;
points=100; %Note to self: the pointval*100 is the problem for different numbers of points
sgrid=linspace(a,b,points);

%[sprobs] = DistFn(['OnePoint'],sgrid,a,b);
[sprobs] = DistFn2('TwoPoints',sgrid,a,b);

%Other parameters for forward solution 
rho=1;
k=1.5;
y0=.2; %note c is normalized so needs to start between 0 and 1
tfinal=20;
tpoints=25;
tspan=linspace(0,tfinal,tpoints); %this is so we can initialize matrix with known number of rows 

%[t, cmat] = ForwardFunction(sgrid, sprobs, rho, k, y0, tspan)
[t, cmat,weightedsol] = RK4FunctionC(sgrid, sprobs, rho, k, y0, tspan);

figure
for i=1:points
    plot(t,cmat(:,i))
    hold on
end
hold on
plot(t,weightedsol,'*')

    figure
    for i=1:length(sgrid)
        hold on
        stem(sgrid(i), sprobs(i),'--o','LineWidth',3,'MarkerSize',12)
    end
    set(gca,"FontSize",20)
    %title('Proportions of Different s Values')
    xlabel('Sensitivity to Treatment (s)')
    ylabel('Proportion of Population')
    %ylim([0,1])

% %plot one curve
% figure
% plot(t,cmat(1,:))

    % %plot all curves
    % figure
    % for i=1:length(sgrid)
    %     hold on
    %     plot(t, cmat(:,i),'-','LineWidth',2)
    %     Legend{i}=strcat('s=',' ', num2str(sgrid(i)));
    % end
    % set(gca,"FontSize",20)
    % %title('Proportions of Different s Values')
    % xlabel('Time')
    % ylabel('Subpopulation Growth')
    % legend(Legend)
    % ylim([0,1])

    %plot weightedsol
    figure
    plot(t,weightedsol,'k','LineWidth',2)
    xlabel('Time')
    ylabel('Aggregated Volume')
    ylim([0,1])
    set(gca,"FontSize",20)

   

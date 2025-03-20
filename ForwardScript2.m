
a=0;
b=1;
points=11;
sgrid=linspace(a,b,points);

[sprobs] = DistFn2('Bigaussian',sgrid,a,b);

%Other parameters for forward solution 
rho=1;
k=1.5;
y0=.2; %note c is normalized so needs to start between 0 and 1
tfinal=10;
tpoints=25;
tspan=linspace(0,tfinal,tpoints); %this is so we can initialize matrix with known number of rows 

[t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);

figure
for i=1:length(sgrid)
    hold on
    stem(sgrid(i), sprobs(i),'--o','LineWidth',3,'MarkerSize',12)
end
set(gca,"FontSize",20)
xlabel('Sensitivity to Treatment (s)')
ylabel('Proportion of Population')

%plot all c(t,s) curves
figure
for i=1:length(sgrid)
    hold on
    plot(t, cmat(:,i),'-','LineWidth',2)
    Legend{i}=strcat('s=',' ', num2str(sgrid(i)));
end
set(gca,"FontSize",20)
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

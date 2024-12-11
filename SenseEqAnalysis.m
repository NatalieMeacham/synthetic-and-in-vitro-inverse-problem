%The goal of this file is just to plot the equilibrium points of the
%sensitive growth model for varying rho, k, and s, and write stability
%analysis code so that I can come back with specific values later and check
%the stability as necessary. 
sgrid=linspace(0,1,15);
i=12;
%rho=2/3;
%k=1;
rho=0.3;
k=0.45;

tspan=linspace(0,60,90);
y0=0.2;
[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*sgrid(i)*csol, tspan, y0); 
%[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)), tspan, y0);
%%version without death term

posorneg=rho/(rho+k);
if sgrid(i)<posorneg 
    disp('The nonzero eq. will be positive')
elseif sgrid(i)>posorneg 
    disp('The nonzero eq will be negative')
else 
    disp('The nonzero eq is also zero in this case')
end

eq1=0;
eq1vec=ones(length(t),1).*eq1;
eq2=(rho*sgrid(i) - rho + k*sgrid(i))/(rho*sgrid(i)-rho);
eq2vec=ones(length(t),1).*eq2;

figure
plot(t,csol,'LineWidth',2)
hold on
plot(t,eq1vec,'LineWidth',2)
hold on
plot(t,eq2vec,'LineWidth',2)
legend('True Soln.' ,'Zero Eq.', 'Nonzero Eq.')
xlabel('Time')
ylabel('Cell Density')
set(gca,"FontSize",20)

%f'(c) at zero eq:
fprime1=rho - rho*sgrid(i) - k*sgrid(i);

%f'(c) at nonzero eq:
fprime2=rho*sgrid(i) - rho + k*sgrid(i);

%check and state stability of first eq
if fprime1 > 0
    disp('zero eq is unstable')
elseif fprime1<0 
    disp('zero eq is stable')
else
    disp('f prime at zero eq. is 0')
end

%check and state stability of second eq
if fprime2 > 0
    C=string(eq2);
    A=strcat('nonzero eq c=',C,' is unstable');
    disp(A)
elseif fprime2<0 
    %disp('nonzero eq is stable')
    F=string(eq2);
    B=strcat('nonzero eq c=',F,' is stable');
    disp(B)
else
    disp('unclear whether second nonzero eq. stable or unstable')
end


points=30;
sgrid=linspace(0,1,points);
eq2fn= @(i) (rho.*sgrid(i) + k.*sgrid(i) - rho)/(rho.*sgrid(i) - rho); 
vec=[];
for i=1:1:points
    vec(i)=eq2fn(i);
end
vec;
figure
plot(sgrid,vec,'LineWidth',2)
hold on
%plot(sgrid,rho/(rho+k),'*','LineWidth',2)
%hold on
%plot(rho/(rho+k),rho/(rho+k),'o','MarkerSize',20)
%hold on
xline(rho/(rho+k),'LineWidth',2)
hold on
plot(sgrid,zeros(1,points),'LineWidth',2)
legend('Nonzero Eq.','Bifurcation','Zero','Location','southwest')
%ylim([-2 2])
axis([0 1 -1 1])
xlabel('Sensitivity to Treatment {\it s}')
ylabel('Nonzero Eq. Point as Function of {\it s}')
%title('Change in Nonzero Eq. as s Changes')
set(gca,"FontSize",20)


%rho=1/4;
%k=3/4;
bifurcationindex = posorneg*100 + 1;
sgrid=linspace(0,1,101);
sgrida = sgrid(1:bifurcationindex);
sgridb=sgrid(bifurcationindex:end);
eq2fn= @(i) (rho.*sgrid(i) + k.*sgrid(i) - rho)/(rho.*sgrid(i) - rho); 
beforebivec=[];
for i=1:101
    beforebivec(i)=eq2fn(i);
end

afterbivec=[];
for i=1:101
    afterbivec(i)=eq2fn(i);
end



vec1=linspace(0,1,101);
vec2=linspace(posorneg,1,101);
figure
plot(vec1(1:bifurcationindex),zeros(1,length(vec1(1:bifurcationindex))),'-b','LineWidth',2) %zero eq before bifurcation
hold on
plot(vec1(1:bifurcationindex),beforebivec(1:bifurcationindex),'--g','LineWidth',2) %nonzero eq before bifurcation
hold on
plot(vec2,zeros(1,length(vec2)),'b--','LineWidth',2) %zero eq after bifurcation
hold on
plot(vec1(bifurcationindex:end),afterbivec(bifurcationindex:end),'g','LineWidth',2) %nonzeroeq after bifurcation
hold on
xline(rho/(rho+k),'LineWidth',2)
% hold on
% plot(sgrid,vec)
axis([0 1 -1 1])
xlabel('Sensitivity to Treatment{\it s}')
ylabel('Nonzero Eq. Point as Function of{\it s}')
%title('Change in Nonzero Eq. as s Changes')
set(gca,"FontSize",20)
legend('Zero Eq. (Unstable)', 'Nonzero Eq. (Stable)', 'Zero Eq. (Stable)', 'Nonzero Eq. (Unstable)', 'Bifurcation','FontSize',16)

% 
% figure
% plot(vec1,zeros(1,length(vec1)),'--b','LineWidth',2)
% hold on
% plot(vec1,vec(1:26),'g','LineWidth',2)
% hold on
% plot(vec2,zeros(1,length(vec2)),'b','LineWidth',2)
% hold on
% plot(vec2,vec(26:101),'--g','LineWidth',2)
% hold on
% xline(rho/(rho+k),'LineWidth',2)
% % hold on
% % plot(sgrid,vec)
% axis([0 1 -1 1])
% xlabel('Sensitivity to Treatment (s)')
% ylabel('Nonzero Eq. Point as Function of s')
% %title('Change in Nonzero Eq. as s Changes')
% set(gca,"FontSize",20)
% legend('Zero Eq. (Unstable)', 'Nonzero Eq. (Stable)', 'Zero Eq. (Stable)', 'Nonzero Eq. (Unstable)', 'Bifurcation')
% 
% %This needs to be updated with new forward function and new equilibria
% 
% a=0;
% b=1;
% points=150; %this used to be 50
% disttype = 'Bigaussian';
% rho=1; %make these inputs 
% k=1.5; %maximal death rate due to treatment %bigger than rho
% y0=.2; %note c is normalized so needs to start between 0 and 1
% tfinal=10;
% tpoints=25;
% tspan=linspace(0,tfinal,tpoints);
% 
% sgrid=linspace(a,b,points);
% [sprobs] = DistFn('Bigaussian',sgrid,a,b);
% 
% [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);
% eq2fn= @(i) (rho.*sgrid(i) + k.*sgrid(i) - rho)/(rho.*sgrid(i) - rho); 
% 
% 
% 
% 
% s=0.3;
% dcdt = @(c) rho.*c.*(1-c).*(1-s) - k.*s.*c;
% cvec = linspace(-1,1,100);
% figure
% plot(cvec,dcdt(cvec),'b','LineWidth',2)
% hold on
% xline(0,'r-','LineWidth',2)
% hold on
% yline(0,'g','LineWidth',2)
% hold on
% xline((rho*s + k*s - rho)/(rho*s - rho),'m-','LineWidth',2)
% xlabel('c(t)')
% ylabel('$\frac{dc}{dt}$','interpreter','latex')
% title('Equilibrium Stability for $s=0.3<\frac{\rho}{\rho + k}$','interpreter','latex')
% legend('$\frac{dc}{dt}$','Zero Eq.','$\frac{dc}{dt}=0$','Nonzero Eq.','interpreter','latex')
% 
% s=0.5;
% dcdt = @(c) rho.*c.*(1-c).*(1-s) - k.*s.*c;
% figure
% plot(cvec,dcdt(cvec),'b','LineWidth',2)
% hold on
% xline(0,'r-','LineWidth',2)
% hold on
% yline(0,'g','LineWidth',2)
% hold on
% xline((rho*s + k*s - rho)/(rho*s - rho),'m-','LineWidth',2)
% xlabel('c(t)')
% ylabel('$\frac{dc}{dt}$','interpreter','latex')
% title('Equilibrium Stability for $s=0.5>\frac{\rho}{\rho + k}$','interpreter','latex')
% legend('$\frac{dc}{dt}$','Zero Eq.','$\frac{dc}{dt}=0$','Nonzero Eq.','interpreter','latex')
% 
% s=rho/(rho + k);
% dcdt = @(c) rho.*c.*(1-c).*(1-s) - k.*s.*c;
% figure
% plot(cvec,dcdt(cvec),'b','LineWidth',2)
% %hold on
% %xline(0,'r-','LineWidth',2)
% hold on
% yline(0,'g','LineWidth',2)
% hold on
% xline((rho*s + k*s - rho)/(rho*s - rho),'m-','LineWidth',2)
% xlabel('c(t)')
% ylabel('$\frac{dc}{dt}$','interpreter','latex')
% title('Equilibrium Stability for $s=0.4=\frac{\rho}{\rho + k}$','interpreter','latex')
% legend('$\frac{dc}{dt}$','$\frac{dc}{dt}=0$','Nonzero Eq.=Zero Eq.','interpreter','latex')

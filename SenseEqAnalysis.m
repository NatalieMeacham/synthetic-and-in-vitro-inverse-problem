%Plot the equilibria of the sensitive growth model for varying 
%rho, k, and s
%Assess stability of the two equilibria  
%Generate supplementary figure 8

sgrid=linspace(0,1,15);
i=14;
rho=0.3; 
k=0.45;

tspan=linspace(0,60,90);
y0=0.2;
[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*sgrid(i)*csol, tspan, y0); 

%check relationship between s and rho/(rho+k)
posorneg=rho/(rho+k);
if sgrid(i)<posorneg 
    disp('The nonzero eq. will be positive')
elseif sgrid(i)>posorneg 
    disp('The nonzero eq will be negative')
else 
    disp('The nonzero eq is also zero in this case')
end

%set up vectors for first plot
eq1=0;
eq1vec=ones(length(t),1).*eq1;
eq2=(rho*sgrid(i) - rho + k*sgrid(i))/(rho*sgrid(i)-rho);
eq2vec=ones(length(t),1).*eq2;

%f'(c) at zero eq:
fprime1=rho - rho*sgrid(i) - k*sgrid(i);

%f'(c) at nonzero eq:
fprime2=rho*sgrid(i) - rho + k*sgrid(i);

%check and state stability of first eq
if fprime1 > 0
    disp('Zero eq. is unstable.')
elseif fprime1<0 
    disp('Zero eq. is stable.')
else
    disp('Value of f prime at zero eq. is 0.')
end

%check and state stability of second eq
if fprime2 > 0
    C=string(eq2);
    A=strcat('Nonzero eq. c=',C,' is unstable.');
    disp(A)
elseif fprime2<0 
    %disp('nonzero eq is stable')
    F=string(eq2);
    B=strcat('Nonzero eq. c=',F,' is stable.');
    disp(B)
else
    disp('Unclear whether second nonzero eq. stable or unstable.')
end

%plot true solution and equilibria
figure
plot(t,csol,'r','LineWidth',2)
hold on
if fprime1>0
    plot(t,eq1vec,'--b','LineWidth',2)
elseif fprime1<0 
    plot(t,eq1vec,'b','LineWidth',2)
else
    disp('Value of f prime at zero eq. is 0.')
end
hold on
if fprime2 > 0
    plot(t,eq2vec,'--g','LineWidth',2)
elseif fprime2<0 
    plot(t,eq2vec,'g','LineWidth',2)
else
    disp('Unclear whether second nonzero eq. stable or unstable.')
end
legend('True Soln.' ,'Zero Eq.', 'Nonzero Eq.')
xlabel('Time')
ylabel('Cell Density')
titlestring=strcat('True Solution for {\it s}=',string(sgrid(i)))
title(titlestring)
set(gca,"FontSize",20)

%set up vectors for bifurcation plot
points=30;
sgrid=linspace(0,1,points);
eq2fn= @(i) (rho.*sgrid(i) + k.*sgrid(i) - rho)/(rho.*sgrid(i) - rho); 
vec=[];
for i=1:1:points
    vec(i)=eq2fn(i);
end

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

%plot bifurcation 
vec1=linspace(0,1,101);
vec2=linspace(posorneg,1,101);
figure
plot(vec1(1:bifurcationindex),zeros(1,length(vec1(1:bifurcationindex))),'--b','LineWidth',2) %zero eq before bifurcation
hold on
plot(vec1(1:bifurcationindex),beforebivec(1:bifurcationindex),'-g','LineWidth',2) %nonzero eq before bifurcation
hold on
plot(vec2,zeros(1,length(vec2)),'b-','LineWidth',2) %zero eq after bifurcation
hold on
plot(vec1(bifurcationindex:end),afterbivec(bifurcationindex:end),'--g','LineWidth',2) %nonzeroeq after bifurcation
hold on
xline(rho/(rho+k),'LineWidth',2)
axis([0 1 -1 1])
xlabel('Sensitivity to Treatment{\it s}')
ylabel('Nonzero Eq. Point as Function of{\it s}')
set(gca,"FontSize",20)
legend('Zero Eq. (Unstable)', 'Nonzero Eq. (Stable)', 'Zero Eq. (Stable)', 'Nonzero Eq. (Unstable)', 'Bifurcation','FontSize',16)


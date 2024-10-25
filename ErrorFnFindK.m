function[err]=ErrorFnFindK(k, data,tspan,IC,weights)

%bestapprox=(gencmatC*s)'; %bestapprox gets weightedsol by multiplying s*csol
%[~,c]=ode45(@(t,c) -k*c,tspan,IC ); %regular
%[~,c]=ode45(@(t,c) -k*c*(1 - c),tspan,IC ); %logistic


rhomax=0.0682;
[~,c]=ode45(@(t,c) -k*rhomax*c*(1 - c),tspan,IC );

%[c]=ForwardFunctionFindRho(rho, data,tspan,IC)

err=sum((weights).*((c' - data)).^2);
    %sum((weights).*((bestapprox - data)).^2);
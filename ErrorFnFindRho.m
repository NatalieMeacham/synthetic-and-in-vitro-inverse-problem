function[err]=ErrorFnFindRho(rho, data,tspan,IC,weights)

[~,c]=ode45(@(t,c) rho*c*(1-c),tspan,IC );

err=sum((weights).*((c' - data)).^2);

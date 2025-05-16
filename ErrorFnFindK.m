function[err]=ErrorFnFindK(k, data,tspan,IC,weights)

[~,c]=ode45(@(t,c) -k*c,tspan,IC );

err=sum((weights).*((c' - data)).^2);
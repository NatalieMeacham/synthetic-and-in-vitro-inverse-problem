function[t,c]=ForwardFnFindK(k,tspan,IC)
rhoval=0.0692;
%[t,c]=ode45(@(t,c) -k*c,tspan,IC ); %regular death term
%[t,c]=ode45(@(t,c) k*c,tspan,IC );
%[t,c]=ode45(@(t,c) -k*c*rhoval,tspan,IC ); %death term with max rho value
[t,c]=ode45(@(t,c) -k*c*(1-c),tspan,IC ); %logistic death term 
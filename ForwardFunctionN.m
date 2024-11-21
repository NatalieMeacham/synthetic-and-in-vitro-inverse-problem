function [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan)
    %Solve ODE for each point on sgrid:
    
    cmat=zeros(length(tspan),length(sgrid));
    %figure
    for i=1:length(sgrid)
        [t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) -k*sgrid(i)*csol, tspan, y0); %regular
       %[t,csol] = ode45(@(t,csol) - k*sgrid(i)*csol, tspan, y0);  %make "all sensitive" model for reviewer 1
        %[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i))- k*rho*sgrid(i)*csol, tspan, y0); %death term has rho
        %[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i))-k*rho*sgrid(i)*csol*(1 - csol), tspan, y0); %logistic death term
        %fill in matrix for csol values over time and different s values
            cmat(:,i)=csol;
  
    end

     %for some time t, dot product the sprobs with t's row in the cmat matrix
     %weightedsol=sprobs*cmat'; %this gives the total weighted c across time
     %this is waht needs to change -> see ER's code 
     weightedsol=(sprobs./(sum(sprobs)))*cmat';
    
 end 
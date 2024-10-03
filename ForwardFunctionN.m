function [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan)
    %Solve ODE for each point on sgrid:
    
    cmat=zeros(length(tspan),length(sgrid));
    %figure
    for i=1:length(sgrid)
        %[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*sgrid(i)*csol, tspan, y0); 
        [t,csol] = ode45(@(t,csol) - k*sgrid(i)*csol, tspan, y0);  %make "all sensitive" model for reviewer 1
        %fill in matrix for csol values over time and different s values
            cmat(:,i)=csol;
        %disp(sgrid(i));
        %disp(csol);
  
         % hold on
         % plot(t,csol,'LineWidth',2)
         % xlabel('time')
         % ylabel('cell density')
         % title('c Solution for Different s Values')
         % hold off
  
    end

     %for some time t, dot product the sprobs with t's row in the cmat matrix
     %weightedsol=sprobs*cmat'; %this gives the total weighted c across time
     %this is waht needs to change -> see ER's code 
     weightedsol=(sprobs./(sum(sprobs)))*cmat';
    
    % N=length(sgrid);
    %  %Legend=cell(N,1)
    %  Legend=cell(N+1,1)
     % for iter=1:N
     % %Future goal: only print a the legend for a few of the curves
     %    Legend{iter}=strcat('csol for s=',' ', num2str(sgrid(iter)));
     % end
     %legend(Legend, 'FontSize',12)
     % hold on
     % y0func = y0*ones(1, length(t));
     % plot(t,y0func,'*','LineWidth',2)
     % hold off
     % %ICstring=strcat('IC=',' ',num2str(y0));
     % Legend{N+1}=strcat('IC=',' ',num2str(y0))
     % legend(Legend, 'FontSize',14)
     % set(gca,"FontSize",20)
     % %need to update the legend so that IC shows as such 
    
     % %plot s distribution
     % figure
     % p=plot(sgrid, sprobs, '*-')
     % p.LineWidth = 3;
     % set(gca,"FontSize",20)
     % title('Proportions of Different s Values')
     % xlabel('sensitivity values')
     % ylabel('proportion of population')
     % ylim([0,1])
   % 
   % 
     % figure
     % q=plot(t, weightedsol)
     % q.LineWidth = 3;
     % %q(2).LineWidth = 2;
     % hold on
     % 
     % r=plot(t,y0func,'*')
     % r.LineWidth = 2;
     % %hold off
     % %ylim([y0 - .1, y0 + .1]) %use this when printing results from different dists
     % ylim([0,1])
     % title('Weighted c Solution')
     % xlabel('time')
     % ylabel('weighted cell density')
     % legend('weighted c','IC')
     % set(gca,"FontSize",20)
% % 
% %     %Z = 4 * ones(1, length(x));
% %     %y0func = y0*ones(1, length(t));
% %     %figure 
% %     %plot(t,y0func)
% %     
% % %     %cmat
 end 
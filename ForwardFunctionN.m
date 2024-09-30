%This should be modified from ForwardFunctionN to ForwardFunctionC, with C
%for competition. 

%use the kc(1 - c)(1-s) version to generate synthetic data called u, and then go forward with plugging in that u to the DE

function [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan)
    %Solve ODE for each point on sgrid:
   

    cmat=zeros(length(tspan),length(sgrid));
%     figure
    for i=1:length(sgrid)
        [t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*sgrid(i)*csol, tspan, y0); %normal version
        %[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*sgrid(i)*csol*(1 - csol), tspan, y0); %logistic version
        %fill in matrix for csol values over time and different s values
        %[t,csol] = ode45(@(t,csol) rho*csol*(1 - csol)*(1 - sgrid(i)) - k*rho*sgrid(i)*csol, tspan, y0); %version where death term is scaled by growth rate
        %[t,csol] = ode45(@(t,csol)  - k*csol, tspan, y0);   %only sensitive version for reviewer 
        cmat(:,i)=csol;
%         %disp(sgrid(i));
%         %disp(csol);

        % hold on
        % plot(t,csol,'LineWidth',2)
        % xlabel('Time')
        % ylabel('Cell Density')
        % %title('c Solution for Different s Values')
        % hold off

    end

     %for some time t, dot product the sprobs with t's row in the cmat matrix
     weightedsol=sprobs*cmat'; %this gives the total weighted c across time %original
     %weightedsol=((max(sprobs)-min(sprobs))/sum(sprobs))*sprobs; %trying to normalize 
     %weightedsol=(sprobs/sum(sprobs))*cmat'; %another attempt to normalize 
   % N=length(sgrid);
   %  Legend=cell(N,1)
   %  %Legend=cell(N+1,1)
   %  for iter=1:N
   %  %Future goal: only print a the legend for a few of the curves
   %     Legend{iter}=strcat('s=',' ', num2str(sgrid(iter)));
   %  end
   %  legend(Legend, 'FontSize',12)
   %  hold on
   %  %y0func = y0*ones(1, length(t));
   %  %plot(t,y0func,'*','LineWidth',2)
   %  %hold off
   %  %ICstring=strcat('IC=',' ',num2str(y0));
   %  %Legend{N+1}=strcat('IC=',' ',num2str(y0))
   %  legend(Legend, 'FontSize',14)
   %  set(gca,"FontSize",20)
   %  %need to update the legend so that IC shows as such 
   % % 
   %  %plot s distribution
   %  figure
   %  %p=plot(sgrid, sprobs, '*-')
   %  p=plot(sgrid, sprobs, '*')
   %  p.LineWidth = 3;
   %  set(gca,"FontSize",20)
   %  %title('Proportions of Different s Values')
   %  xlabel('Sensitivity to Treatment (s)')
   %  ylabel('Proportion of Population')
   %  ylim([0,1])
   %  ylim([0,0.2])
% 
%    % 
%    % 
    % figure
    % q=plot(t, weightedsol)
    % q.LineWidth = 3;
    % %q(2).LineWidth = 2;
    % hold on

%     r=plot(t,y0func,'*')
%     r.LineWidth = 2;
%     
%     %hold off
%     %ylim([y0 - .1, y0 + .1]) %use this when printing results from different dists
%     ylim([0,1])
%     %title('Weighted c Solution')
%     xlabel('Time')
%     ylabel('Aggregated Tumor Volume')
% %     legend('Aggregated Volume','IC')
%     set(gca,"FontSize",20)
%     hold off
% % % 
% % %     %Z = 4 * ones(1, length(x));
% % %     %y0func = y0*ones(1, length(t));
% % %     %figure 
% % %     %plot(t,y0func)
% % %     
% % % %     %cmat
%  end 
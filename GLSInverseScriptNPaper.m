 function [sgrid,sprobs, rsgrid,optweightfromAIC,t,propdata,sweightedsol,cdfS,cdfR]= GLSInverseScriptNPaper(disttype, tpoints, noisesize)
    
 %Create synthetic data and run the inverse problem on that data.
 %Generate figures showing recovered fit, PMF, and CDF versus the original
 %data, PMF/PDF, and CDF respectively. 

 %INPUTS:
    %disttype: distribution type for initial sensitivity
    %distribution
    %tpoints: number of time points in synthetic data
    %noisesize: proportion up to which forward solution is scaled to get
    %noisy synthetic data 

 %OUTPUTS
    %sgrid: sensitivity mesh, determined by points, a, and b
    %sprobs: sensitivity distribution, determined by sgrid and disttype
    %rsgrid: sensitivity mesh for optimal recovered distribution
    %optweightfromAIC: optimal recovered distribution
    %t: time vector output from forward problem
    %propdata: synthetic data
    %sweightedsol: aggregated tumor volume curve from forward problem using
    %optimal recovered distribution
    %cdfS: CDF from original distribution (sprobs)
    %cdfr: CDF from recovered distribution (optweightfromAIC)

    %declare variables 
    points=101; 
    tfinal=50;
    tspan=linspace(0,tfinal,tpoints);
    pointsstr=string(points);
    noisestr=string(noisesize);

    %parameters for synthetic data
    rho=0.3;
    k=0.45;
    y0=.2;
    
    %define sensitivity mesh
    a=0;
    b=1;
    sgrid=linspace(a,b,points);
    
    %create original distribution, synthetic data, add noise
    [sprobs] = DistFn2(disttype,sgrid,a,b);
 
    %generate synthetic data with proportional error
    [t, cmat,weightedsol] = ForwardFunctionN(sgrid, sprobs, rho, k, y0, tspan);
    propdata=weightedsol.*(1 + noisesize*randn(size(weightedsol)));
    
    %loop through IP over different meshes -> find optimal mesh
    rpointsvec=4:1:30;
    AICvec=zeros(length(rpointsvec),1);
    for i=1:length(rpointsvec)
        [gls_optpar,~,AIC_GLS,~,~,~,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,rpointsvec(i),propdata);
        AICvec(i)=AIC_GLS;
    end
    
    [M,I]=min(AICvec);
    optrpoints=rpointsvec(I);
    
    %get optweight from OLSInverseFn for optrpoints
    [optweightfromAIC,~,~,sgrid,sprobs,weightedsol,rsgrid]=GLSInverseFnN(points,disttype,rho,k,y0,tfinal,tpoints,noisesize,optrpoints,propdata);
    
    figure
    plot(rpointsvec,AICvec,'LineWidth',2,'Color','blue')
    xlabel('Number of Points in Recovered Dist.')
    ylabel('AIC Score')
    set(gca,"FontSize",20)
    hold on
    plot(rpointsvec(I),M,'m*','LineWidth',2,'MarkerSize',20)
    hold off
    legend('AIC','Min. AIC','Location','northeast')
    
    %compare output of fmincon with original distribution
    figure
    yyaxis left
    stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) 
    ylabel('Recovered Proportion of Population')
    limsy=get(gca,'YLim');
    set(gca,'Ylim',[0 limsy(2)]);
    ylimval=max(max(sprobs),max(optweightfromAIC));
    if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
        ylim([0 ylimval])
    end
    if strcmp(disttype,'Uniform') == 1
        ylim([0 2.5*median(optweightfromAIC)])
    end
    hold on
    yyaxis right
    plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
    legend('Recovered','Original','Location','northeast')
    xlabel('Sensitivity to Treatment {\it s}')
    ylabel('Original Proportion of Population')
    if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
        ylim([0 ylimval])
    end
    if strcmp(disttype,'Uniform') == 1
        ylim([0 2.5*median(sprobs)])
    end
    set(gca,"FontSize",20)
    
    %make and plot cdfs to accommodate continuous distributions:
    cdfR=zeros(optrpoints,1);
    cdfS=zeros(points,1);
    if strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Uniform') ==1 || strcmp(disttype,'Bigaussian') ==1
        partialsumS=0;
        cdfS(1)=(1/points)*sprobs(1);
        for i=2:points
            partialsumS=trapz(sgrid(1:i),sprobs(1:i));
            cdfS(i)=partialsumS;
        end
        partialsumR=0;
        for i=1:optrpoints
            partialsumR=sum(optweightfromAIC(1:i));
            cdfR(i)=partialsumR;
        end
        disp('Here is a cdf for a continuous dist.')
    elseif strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') ==1
        partialsumS=0;
        for i=1:points
            partialsumS=sum(sprobs(1:i));
            cdfS(i)=partialsumS;
        end
        partialsumR=0;
        for i=1:optrpoints
            partialsumR=sum(optweightfromAIC(1:i));
            cdfR(i)=partialsumR;
        end
        disp('Here is a cdf for a discrete dist.')
    end
    
    figure
    X=rsgrid;
    Y=cdfR;
    stairs(X,Y,'--o','MarkerSize',8,'LineWidth',2,'Color','blue')
    ylim([0 1])
    hold on
    plot(sgrid, cdfS,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
    xlabel('Sensitivity to Treatment {\it s}')
    ylabel('Cumulative Recovered Proportion')
    legend('Recovered','Original','Location','Southeast')
    set(gca,"FontSize",20)
    
    %plot synthetic data versus recovered fit 
    [t,~,sweightedsol] = ForwardFunctionN(rsgrid, optweightfromAIC', rho, k, y0, tspan);
    figure
    plot(t,sweightedsol,'LineWidth',2,'Color','blue')
    hold on
    plot(t,propdata','rd','LineWidth', 2,'MarkerSize',8)
    hold off
    legend('Est. Tumor Volume','Simulated Data','Location','southeast')
    xlabel('Time')
    ylabel('Aggregated Tumor Volume')
    set(gca,"FontSize",20)
    
    %compare output of fmincon with original distribution on same scaled axes 
    if strcmp(disttype,'Uniform') == 1 || strcmp(disttype,'Normal') == 1 || strcmp(disttype,'Bigaussian') == 1
        figure
        stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) 
        ylabel('Proportion of Population')
        hold on
        plot(sgrid,sprobs*length(sprobs)/(sum(sprobs)*length(optweightfromAIC)),'--o','MarkerSize',8,'LineWidth',2,'Color','red')
        limsy=get(gca,'YLim');
        set(gca,'Ylim',[0 limsy(2)]);
        legend('Recovered','Original','Location','northeast')
        xlabel('Sensitivity to Treatment {\it s}')
        set(gca,"FontSize",20)
    end

    if strcmp(disttype,'OnePoint') == 1 || strcmp(disttype,'TwoPoints') == 1 
        figure
        stem(rsgrid,optweightfromAIC,'--b','MarkerSize',8,'LineWidth',2) 
        limsy=get(gca,'YLim');
        set(gca,'Ylim',[0 limsy(2)]);
        hold on
        plot(sgrid,sprobs,'--o','MarkerSize',8,'LineWidth',2,'Color','red')
        legend('Recovered','Original','Location','northeast')
        xlabel('Sensitivity to Treatment {\it s}')
        ylabel('Proportion of Population')
        ylimval=max(max(sprobs),max(optweightfromAIC));
        ylim([0 ylimval])
        set(gca,"FontSize",20)
    end

end

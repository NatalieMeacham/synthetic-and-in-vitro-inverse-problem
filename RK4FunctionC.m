function [h, k1, k2, k3, k4, cmat,weightedsol,uvec,sprobsmat] = RK4FunctionC(tpoints, tfinal, rho, k, sgrid, sprobs, points, y0)
%cmat=zeros(tpoints,points)
tspan=linspace(0,tfinal,tpoints);
%sprobsthingie?
%uinit=; %initial u
h=tspan(2) - tspan(1); %delta t

dcdt = @(csol,u) rho.*csol.*(1 - u).*(1 - sgrid') - k.*sgrid'.*csol

sprobsmat=zeros(points, tpoints);
sprobsmat(:,1)=sprobs;
uvec=zeros(tpoints,1);
y0vec=(y0*ones(points,1))'.*sprobs;
uvec(1)=sum(y0vec);
%uvec(1)=sum((y0*ones(points,1))'.*sprobs)

%u=sum(csol'.*sprobs)
%sprobs
%sprobs=csol'.*sprobs/u

cmat=zeros(length(sgrid),length(tspan));
%cmat(:,1)=y0*ones(length(sgrid),1); %old version
cmat(:,1)=y0vec; %new scaled version
%for loop over t
for j=2:tpoints
    %define k steps
    %k1 = hf(x0, y0)
    %k2 = hf[x0 + (½)h, y0 + (½)k1]
    % k3 = hf[x0 + (½)h, y0 + (½)k2]
    % k4 = hf(x0 + h, y0 + k3)
%keyboard 
    u=uvec(j-1);

    k1=h.*dcdt(cmat(:,j-1),u);
    k2=h.*dcdt(cmat(:,j-1),u + k1/2);
    k3=h.*dcdt(cmat(:,j-1),u + k2/2);
    k4=h.*dcdt(cmat(:,j-1),u + k3);
    cmat(:,j)=cmat(:,j-1) + (k1 + 2*k2 + 2*k3 + k4)/6;%*ones(length(sgrid),1);
    %y1 = y0 + (⅙) (k1 + 2k2 + 2k3 + k4)

    %update u and sprobs
    %u=sum(cmat(:,j)'.*sprobs); %OLD VERSION
    %sprobs
    %sprobs=cmat(:,j)'.*sprobs/u; %OLD SPROBS
    %sprobs=cmat(:,j)'./u; %TRYING SMTH NEW HERE


    %sprobs=cmat(:,j)'./sum(cmat(:,j)); %TRYING SMTH ELSE NEW
    %update uvec
    %uvec(j)=u;
    uvec(j)=sum(cmat(:,j));
    %update sprobsmat
    %sprobsmat(:,j)=sprobs;
end
weightedsol=1;

end
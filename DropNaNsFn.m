function [meanmatnorm,isnanmat, tspanvec, tspanmat,concvecA]=DropNaNsFn(concvec,meanmatnorm,tspan)
%meanmatnorm(any(isnan(meanmatnorm),2),:) = [];
%replace any NaNs row value with 1s
% TFvec=zeros(length(concvec),1);
% TFmat=zeros(length(concvec),length(tspan));
% for a=concvec
%     for i=1:length(tspan)
%     TFmat(a,i)=anynan(meanmatnorm(a,i))
%     %TFvec=[TFvec, TF]
%     end
% end
% 
% for a=concvec
%     for i=1:length(tspan)
%         if TFmat(a,i)==1
%         meanmatnorm(a,i)=[];
%         %disp('updated meanmatnorm with ones')
%         else 
%         %disp('chill')
%         end
%     end
% end

% for a=concvec
%     for t=1:length(tspan)
%         meanmatnorm(isnan(meanmatnorm(a,t))) = [];
%     end
% end

%new tspanvec
isnanmat=zeros(size(meanmatnorm));
tspanvec=zeros(length(concvec),1);
tspanmat=zeros(length(concvec),length(tspan));

for a=concvec
    tspanmat(a,:)=tspan;
end

for a=concvec
    isnanmat(a,:)=isnan(meanmatnorm(a,:));
    %[m,i]=max(isnanmat(a,:));
    k = find(isnanmat(a,:));
    if isempty(k)==1
        tspanvec(a)=14;
    else
        tfillvec=tspanmat(a,:);
        tfillvec(k)=[];
        %tspanmat(a,k)=[];
        tspanmat(a,1:14-length(k))=tfillvec;
        tspanmat(a,15-length(k):end)=zeros;
        tspanvec(a)=14-length(k);
    end
    %tspanmat(a,k)=[];
    B = rmmissing(meanmatnorm(a,:));
    meanmatnorm(a,1:length(B))=B;
    meanmatnorm(a,(length(B)+1):end)=zeros(1,length(tspan)-length(B));
end

%notes on who is missing what nans:
%R250: 2-8 in row 1 are missing (7/14=50%)
%R500: 2 in row 11 is missing (1/11<50%)
%S500: no nans
%S100: 1-14 in row 11 are missing (100%)

%for data missing too many nans, remove that row and adjust all
%notation/dimensions
for a=concvec
    if sum(isnanmat(a,:))>length(isnanmat(a,:))/2
        meanmatnorm(a,:)=[];
        tspanvec(a)=[];
        tspanmat(a,:)=[];
        concvec(a)=[];
        concvecA=concvec;
    else
        disp('We have enough non-nans')
        concvecA=concvec;
    end
end






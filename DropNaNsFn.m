function [meanmatnorm,isnanmat, tspanvec, tspanmat,concvecA]=DropNaNsFn(concvec,meanmatnorm,tspan)

isnanmat=zeros(size(meanmatnorm));
tspanvec=zeros(length(concvec),1);
tspanmat=zeros(length(concvec),length(tspan));

for a=concvec
    tspanmat(a,:)=tspan;
end

for a=concvec
    isnanmat(a,:)=isnan(meanmatnorm(a,:));
    k = find(isnanmat(a,:));
    if isempty(k)==1
        tspanvec(a)=14;
    else
        tfillvec=tspanmat(a,:);
        tfillvec(k)=[];
        tspanmat(a,1:14-length(k))=tfillvec;
        tspanmat(a,15-length(k):end)=zeros;
        tspanvec(a)=14-length(k);
    end
    B = rmmissing(meanmatnorm(a,:));
    meanmatnorm(a,1:length(B))=B;
    meanmatnorm(a,(length(B)+1):end)=zeros(1,length(tspan)-length(B));
end

%for data missing too many NaNs, remove that row and adjust all dimensions
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






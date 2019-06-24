function Hurst=HurstCompute(Xtimes)
%HurstCompute

LengthX=length(Xtimes);

[FactorMatrix,FactorNum]=HurstFactorization(LengthX);

LogRS=zeros(FactorNum,1);

LogN=zeros(FactorNum,1);

for i=1:FactorNum
    
    dataM=reshape(Xtimes,FactorMatrix(i,:));
    MeanM=mean(dataM);
    SubM =dataM-repmat( MeanM,FactorMatrix(i,1),1) ;
    RVector=zeros(FactorMatrix(i,2),1);
    SVector=zeros(FactorMatrix(i,2),1);
 
    for j=1:FactorMatrix(i,2)
        SubVector=cumsum(SubM(:,j));
        RVector(j)=max(SubVector)-min(SubVector);
        SVector(j)=std(dataM(:,j));
    end
    LogRS(i)=log(sum( RVector./SVector)/ FactorMatrix(i,2));
    LogN(i)=log(FactorMatrix(i,1));
end

HurstExponent=polyfit(LogN,LogRS,1);
Hurst = HurstExponent(1);

function  [FactorMatrix,FactorNum]=HurstFactorization(x)
%hurstFactorization
N=floor(x/4);
FactorNum=0;
for i=4:N
    if mod(x,i)==0
       FactorNum=FactorNum+1;
       FactorMatrix(FactorNum,:)=[i,x/i];
    end
end
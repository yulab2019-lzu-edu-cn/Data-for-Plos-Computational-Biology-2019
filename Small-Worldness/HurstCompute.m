function HurstExponent = HurstCompute(Xtimes)
%HurstCompute
%�������ΪXtimes
LengthX = length(Xtimes);
%������ʽ�ֽ�
[FactorMatrix, FactorNum] = HurstFactorization(LengthX);
%����LogRS��Ϊ�����������ĳ�ʼһ��Ϊ0
LogRS = zeros(FactorNum, 1);
%����LogN
LogN = zeros(FactorNum, 1);
%�������
for i = 1 : FactorNum
    %������ʽ�ֽⷽ�������������з���
    %���� FactorMatrix(i,:)=[8    30]
    %��240��Ԫ�ص���������ת��Ϊ8X30�ľ���
    dataM = reshape(Xtimes, FactorMatrix(i, :));
    %�������ÿ�еľ�ֵ
    MeanM = mean(dataM);
    %ִ��
    SubM = dataM - repmat(MeanM, FactorMatrix(i, 1), 1) ;
    RVector = zeros(FactorMatrix(i, 2), 1);
    SVector = zeros(FactorMatrix(i, 2), 1);
   %���㣨R/S��n���ۼ�
    for j = 1 : FactorMatrix(i, 2)
        %SubVector=zeros(FactorMatrix(i,1),1);
        SubVector = cumsum(SubM(:, j));
        RVector(j) = max(SubVector) - min(SubVector);
        SVector(j) = std(dataM(:, j));
    end
    %�ֱ����LogRS��LogN
    LogRS(i)=log( sum( RVector./SVector)/ FactorMatrix(i,2) );
    LogN(i)=log( FactorMatrix(i,1) );
end
%ʹ����С���˷����лع飬�����˹��ָ��HurstExponent
HurstExponent=polyfit(LogN,LogRS,1);

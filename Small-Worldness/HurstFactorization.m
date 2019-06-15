function  [FactorMatrix, FactorNum] = HurstFactorization(x)
%hurstFactorization
%���ӷֽ�, ��4��ʼ��X/4����
%floor������ʾ��������
N = floor(x / 4);
%����������ʼΪ0
FactorNum = 0;
%���ӷֽ�, ��4��ʼ��X/4����
for i = 4 : N
    %i���Ա�x����,���õ�һ��ֽⷽ��
    if mod(x, i) == 0
       %��������+1
        FactorNum = FactorNum+1;
        %�����з����洢��FactorMatrix��
        FactorMatrix(FactorNum, :) = [i, x / i];
    end
end
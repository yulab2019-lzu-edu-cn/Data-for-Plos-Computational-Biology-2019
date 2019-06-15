function  [FactorMatrix, FactorNum] = HurstFactorization(x)
%hurstFactorization
%因子分解, 以4开始以X/4结束
%floor函数表示四舍五入
N = floor(x / 4);
%方案数量初始为0
FactorNum = 0;
%因子分解, 以4开始以X/4结束
for i = 4 : N
    %i可以被x整除,即得到一组分解方案
    if mod(x, i) == 0
       %方案数量+1
        FactorNum = FactorNum+1;
        %将可行方案存储到FactorMatrix中
        FactorMatrix(FactorNum, :) = [i, x / i];
    end
end
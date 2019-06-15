function [Ge, Ge_r, Cost, Cost_r, sigma] = netpropv2_td(eij0, td)%�������ϵ���������ֵ
%��ز�������� Disrupted small-world networks in schizophrenia
eij0 = abs(eij0);
n = size(eij0, 1);
eij = eij0;
eij(eij <= td) = 0;          %ȡ��ֵ
eij(eij > td) = 1;

for i = 1:n
    eij(i, i) = 0;
end

srand = sym_generate_srand(eij);



P = pmin(eij);                   %�����������·��
P_r = pmin(srand);               %�����������·��_�������
ee = 1./P;
ee_r = 1./P_r;

for i = 1:n
    ee(i, i) = 0;
    ee_r(i, i) = 0;
end

nl = n;
Lp = sum(sum(P))/(nl.*(nl - 1));  %ƽ�����·��
Lp_r = sum(sum(P_r))/(nl.*(nl - 1));  %ƽ�����·��_�������
lambda = Lp / Lp_r;

for i = 1 : n                      %����ÿ���ڵ��������
    gi = eij;
    A = find(gi(i, :) == 0); 
    gi(:, A) = [];
    gi(A, :) = [];   
    Pi = pmin(gi);
    ni = size(Pi, 1);
    
    gi_r = srand;
    A_r = find(gi_r(i, :) == 0); 
    gi_r(:, A_r) = [];
    gi_r(A_r, :) = [];   
    Pi_r = pmin(gi_r);
    ni_r = size(Pi_r, 1);
    
    Ci(1, i) = (sum(sum(gi)))./(ni.*(ni - 1));%ÿ��������ľ���ϵ��
    Ci_r(1, i) = (sum(sum(gi_r)))./(ni_r.*(ni_r - 1));%ÿ��������ľ���ϵ��_�������
end

Cp = mean(Ci);
Cp_r = mean(Ci_r);
gamma = Cp / Cp_r;
sigma = gamma / lambda;

Ge = sum(sum(ee)) /(n * (n - 1));    %global efficiency
Ge_r = sum(sum(ee_r)) /(n * (n - 1));    %global_r efficiency

Cost = sum(sum(eij)) / (n * (n - 1)); %Kcost
Cost_r = sum(sum(srand)) / (n * (n - 1)); %Kcost_r
end


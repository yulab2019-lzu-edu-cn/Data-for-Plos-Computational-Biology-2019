clc;clear all

load DTI_gunter
%Branching parameters at different model thersholds
%tm            GH model threshold
%theta         Branching parameter
%DTI_gunter    DTI structural connection matrix     

tm = [0.01,0.1,0.2,0.25,0.3,0.35,0.4,0.45,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.65,0.7,0.75,0.8,0.9,1];
nn = 100;  %simulation times


for i = 1:length(tm)
    ratio = zeros(1,nn);
    tt = tm(i);
    for k = 1:nn
        
        D=GHmodel_mex(tt,0.995,0.98,55,DTI_gunter); %time series of the GH model
        D(D<=0)=0;
        D(1:10000,:) = []; %delete the first 10000 steps for stable activity
       
        %steps were binned to accord with the result from avalanche distribution, the bin width here is 2 steps 
        S0 = sum(D,2);
        S = zeros(14000,1);
        for tt0 = 1:14000
            ts = (tt0-1).*2+1;
            S(tt0,1) = sum(S0(ts:(ts+1),:));
        end
        
        
        %Estimating branching parameter
        
        %delete continuous steps with no activity
        index0=find(S==0);
        ind_diff=index0(2:length(index0))-index0(1:length(index0)-1);
        ko=find(ind_diff==1)+1;
        S(index0(ko))=[];
        
        %branching in each avalanche
        n = length(S);
        t1 = 0;
        for j = 2:n-1
            if S(j-1)==0&&S(j)~=0       %estimated only from the first and second timebin of an avalanche.
                t1 = t1+1;
                d(t1,1) = S(j+1)./S(j);         
            end
        end
        ratio(1,k) = sum(d)./t1;
        
    end
    theta(1,i) = mean(ratio);
end

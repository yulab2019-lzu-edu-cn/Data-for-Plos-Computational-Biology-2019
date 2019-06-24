%Caculating the hurst exponent and network properties of experimental data

datanumber = linspace(490001,490100);
td = 0.4;              %network connection threshold
aa = [3 25 36 50 94];  
datanumber(aa) = [];   %subjects with severe head movement were removed from analysis

for i = 1:95
    
    datan = datanumber(i);
    FileName = ['ROISignals_' num2str(datan)];
    load (FileName);
    data = ROISignals(:,1:90);
    
    corrmatrix = corr(data);
    
    % fisher z-transformation
    corrmatrix = log((1+corrmatrix)./(1-corrmatrix))/2;
    for l = 1:90
        corrmatrix(l,l) = 1;
    end
    
    [output(1,i),output(2,i),output(3,i),output(4,i),output(5,i),output(6,i),output(7,i)] = netprop(corrmatrix,td);
    
%    Hurst exponents were calculated for each mean BOLD time series
     meandata = mean(data,2);   
     output(8,i) = HurstCompute(meandata);

end

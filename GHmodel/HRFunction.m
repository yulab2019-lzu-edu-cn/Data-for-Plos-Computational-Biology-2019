function  HRF=HRFunction()
% where o is the onset delay, d is the time-scaling, and p is an integer phase-delay (the peak delay is given by
% pd, and the dispersion by pd 2  This function is bounded and pos- itively skewed (unlike a Gaussian for example).
g=0:0.01:10;
o=0;  a=0.6; p=3;
HRF=((g-o)./a).^(p-1).*((exp(-(g-o)./a))./a.*factorial(p-1)).^(p-1);

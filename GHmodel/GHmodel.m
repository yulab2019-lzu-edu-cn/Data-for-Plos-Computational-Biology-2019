function [D]=GHmodel(T,rr1,rr2,delay,DTI_gunter)

% DTI_gunter           the  DTI structural connection matrix     
% T                    excitation threshold
% rr1=0.005            self-excitation rate
% rr2=0.02;            R¡úQ rate after the delayed steps
% delay=55;            delayed steps, the shortest refractory state length
% D                    State of brain regions at each step:
%                      1  Excited state
%                      0  Quite state
%                      <0 Refractory period&state


DTI_gunter=DTI_gunter*0.01;

D=zeros(300001,90);

D(1,:)=uint8(rand(1,90));              %intiating the model with random state

for ll=1:300000                             
    for k=1:90                         %90 brain regions
        switch D(ll,k)
            case 1
                D(ll+1,k)=-delay;
            case 0
                r_1=rand(1);
                Sum=0;                
                for j=1:90
                    if D(ll,j)>=0
                        Sum=Sum+D(ll,j)*DTI_gunter(k,j);    
                    else
                        continue
                    end
                end
                
                if Sum>=T             %reaching the excitation & model threshold Tm
                    D(ll+1,k)=1;
                elseif r_1>=rr1
                    D(ll+1,k)=1;
                else
                    D(ll+1,k)=0;
                end
                
            case -1  
                r_2=rand(1);
                
                if r_2>=rr2
                    D(ll+1,k)=0;
                else
                    D(ll+1,k)=-1;
                end
            otherwise
                D(ll+1,k)=D(ll,k)+1;
        end
    end
end

end


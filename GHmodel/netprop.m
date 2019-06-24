function [ge,le,L,gc,ce,S,nn] = netprop(eij0,td)
%Generate a network and caculat its topological properties from a correlation matrix
% eij0  correlation matrix
% td    connection threshold

% ge    Global efficiency 
% le    Local efficiency
% L     Characteristic path length
% gc    Clustering coefficient
% ce    Connection strength
% S     Sparsity
% nn    Number of nodes not connected to the network

eij0 = abs(eij0);
eij = eij0;
eij(eij<= td) = 0;          
eij(eij>td) = 1;
P = pmin(eij);                   %shortest path between nodes in the network

%remove nodes not connected to the network
a = zeros(1,90);
for i = 1:90
a(i) = length(find(P(i,:) == Inf));
end
b = find(a>45);
eij0(b,:) = [];
eij0(:,b) = [];
nn = length(b);

n=size(eij0,1);
eij = eij0;
eij(eij<= td) = 0;          
eij(eij>td) = 1;
P = pmin(eij);                   
ee = 1./P;

for i = 1:n
    ee(i,i) = 0;
end

ge = sum(sum(ee))./(n*(n-1));    %global efficiency

nl = n;
L = sum(sum(P))/(nl.*(nl-1));  

gc_i = zeros(1,n);
le_i = zeros(1,n);
ce_i = zeros(1,n);

for i = 1:n
    eij(i,i) = 0;
end

for i = 1:n                      %sub network of each node
    gi = eij;
    A = find( gi(i,:) == 0 ); 
    gi(:,A) = [];
    gi(A,:) = [];   
    Pi = pmin(gi);
    ni = size(Pi,1);
    
    ce_i(i) = sum((eij(i,:).*eij0(i,:)))/ni;
    
    gc_i(1,i) = (sum(sum(gi)))./(ni.*(ni-1));
    ee = 1./Pi;
    for j = 1:ni
        ee(j,j) = 0;
    end
    le_i(1,i) = sum(sum(ee))./(ni*(ni-1));%local efficiency
end
le_i(isnan(le_i)) = 0;
le = mean(le_i);
gc_i(isnan(gc_i)) = 0;  
gc = mean(gc_i);
ce_i(isnan(ce_i)) = 0;
ce = mean(ce_i);
S = sum(sum(eij))/(n.*(n-1)); 

end

function [P] = pmin(P)
%shortest path length

n=size(P,1);
P(P==0) = Inf;
P = P-diag(diag(P));

for k=1:n
   for i=1:n
      for j=1:n
         if P(i,k)+P(k,j)<P(i,j)
            P(i,j)=P(i,k)+P(k,j);
         end
      end
   end
end
for i = 1:n
    P(i,i) = 0;
end
end

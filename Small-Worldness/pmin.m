function [P] = pmin(P)%最短路径算法 P是走过的连接数 path是路径矩阵，暂时不用
n=size(P, 1);
P(P == 0) = Inf;
P = P - diag(diag(P));
% path=zeros(n,n);
% for i = 1 : n
%    for j = 1 : n
%       if P(i , j)~=inf
%          path(i , j) = j;
%       end
%    end
% end
for k = 1 : n
   for i = 1 : n
      for j = 1 : n
         if P(i, k) + P(k, j) < P(i, j)
            P(i, j) = P(i, k) + P(k, j);
%             path(i, j) = path(i, k);
         end
      end
   end
end
for i = 1 : n
    P(i, i) = 0;
end

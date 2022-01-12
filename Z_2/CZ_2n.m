function f=CZ_2n(g,c)
%f=sum(c(i)*Qindex2Z2(i,f.n))=sum(c(i)*g(i,:)), c(i) is complex number ,
% let group be Z_2^n, then g \in N^{m * n} means the term [x_1^g(1,1)*x_2^g(1,2)*...x_n^g(1,n);...;x_1^g(m,1)*x_2^g(m,2)*...x_n^g(m,n)]
if nargin==1
    f.c=sparse(2^g,1);
    f.n=g;
else
    n=length(g(1,:));m=size(g,1);
    f.n=n;
    Index=zeros(m,1);
    Indexy=ones(m,1);
    for i=1:m
        Index(i)=Z22Qindex(g(i,:));
    end
    N=2^n;
    f.c=sparse(Index,Indexy,c,N,1);
end
end

function k=Z22Qindex(g)
% {1,2}^n \mapsto N
g=mod(g,2);
g(g==0)=2;
g=(g(1)<47)*47+g;
g=char(g);
g=g(end:-1:1);
k=bin2dec(g);
k=k+1;
end
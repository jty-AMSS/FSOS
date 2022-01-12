function f=CZprod(f1,f2)
% prod the product of f1 and f2

n=f1.n;
N=2^n;
A=f1.c;g1=find(A);
B=f2.c;g2=find(B);
C=zeros(N,1);
for i=1:length(g1)
    for j=1:length(g2)
        gt=Qindex2Z2(g1(i),n);
        ht=Qindex2Z2(g2(j),n);
        g=mod(gt+ht,2);
        g(g==0)=2;
        g=Z22Qindex(g);
        C(g)=C(g)+A(g1(i))*B(g2(j));
    end
end
f.c=sparse(C);
f.n=n;
end
function f=CZsum(f1,f2)
% prod the product of f1 and f2
n=f1.n;
f.c=sparse(f1.c+f2.c);
f.n=n;
end
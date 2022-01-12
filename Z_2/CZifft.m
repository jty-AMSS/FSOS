function V=CZifft(f)
%ifftn: C[Z_2^n]->C^{2*2*2*2..*2 }( n times 2) s.t.
%V(x0,x1,...,xn)=f(x1,x2,...,xn) 
% here x0 is the element in add -group that is 
%f(-1,1,-1,-1)=v(2,1,2,2)
g=find(f.c);
N=f.n;
N=ones(1,f.n)+1;
T=zeros(N);
A=f.c;
for i=1:length(g)
    t=g(i);
    t=Qindex2Z2(t,f.n);
    t=3-t;
    t=Z22Qindex(t);
    T(t)=A(g(i));
end
v=ifftn(T);
for i=1:2^(f.n)
    t=Qindex2Z2(i,f.n);
    t=3-t;
    t=Z22Qindex(t);
    V(t)=v(i);
end
V=2^(f.n)*reshape(V,N);
end
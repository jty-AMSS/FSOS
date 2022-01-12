function h=CZfft(v)
%ifftn: C[Z_2^n]->C^{2*2*2*2..*2 }( n times 2) s.t.
%V(x0,x1,...,xn)=f(x1,x2,...,xn) 
% here x0 is the element in add -group that is 
%f(-1,1,-1,-1)=v(2,1,2,2)
g=find(v);
N=size(v);
n=length(N);
T=zeros(N);
v=v/(numel(v))^2;
A=v;
for i=1:length(g)
    t=g(i);
    t=Qindex2Z2(t,n);
    t=3-t;
    t=Z22Qindex(t);
    T(t)=A(g(i));
end
v=fftn(T);
for i=1:2^(n)
    t=Qindex2Z2(i,n);
    t=3-t;
    t=Z22Qindex(t);
    f(t)=v(i);
end
f=2^(n)*reshape(f,N);
h.c=sparse(f(:));
h.n=n;
end
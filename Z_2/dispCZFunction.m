function g=dispCZFunction(f)
A=f.c;n=f.n;N=2^n;
x=find(A);
g=@(y)0;
for i=1:length(x)
    t=x(i);
    t=Qindex2Z2(t,n);
    t(t==1)=1;
    t(t==2)=0;
    g=@(y)g(y)+A(x(i))*prod(y(:).^t(:));
end
end
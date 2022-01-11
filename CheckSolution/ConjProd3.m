function c=ConjProd3(A,B)
N=size(A);
c=zeros(N);
[x,y,z]=find(A);
[q,w,e]=find(B);
for i=1:length(x)
    for j=1:length(q)
        a=mod(q(j)-x(i),N(1));
        a(a==0)=N(1);
        b=mod(w(j)-y(i),N(2));
        b(b==0)=N(2);
        c(a,b)=c(a,b)+conj(z(i))*e(j);
    end
end
end
function Qindex=GenerateQindex(c)
N=size(c);
Qindex=zeros(prod(N),length(N));
%Q=sym(Q);×¢ÊÍµô¼ÓËÙ
for i=1:prod(N)
    K=LexOrderRank2vector(i,N);
    Qindex(i,:)=K+1;
end
end


function i=LexOrderRank2vector(z,N)
n=length(N);
z=z-1;
for k=n:-1:1
    i(k)=mod(z,N(k));
    z=(z-i(k))/N(k);
end
end
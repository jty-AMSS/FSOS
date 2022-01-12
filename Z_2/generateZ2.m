function X=generateZ2(n,type)
%输出 ： 2^n * n 的矩阵 ，每一行为一个Z_2^n的元素 默认为{-1,1}^n, type>=0的时候是 {0,1}^n
if nargin==1
    type=-1;
end
if type<0
    X=zeros(2^n,n);
    for i=0:size(X,1)-1
        t=dec2bin(i,n);
        t=2*(double(t)-48.5);
        X(i+1,:)=t;
    end
else
    X=zeros(2^n,n);
    for i=0:size(X,1)-1
        t=dec2bin(i,n);
        t=(double(t)-48);
        X(i+1,:)=t;
    end
end

end
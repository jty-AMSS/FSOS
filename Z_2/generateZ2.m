function X=generateZ2(n,type)
%��� �� 2^n * n �ľ��� ��ÿһ��Ϊһ��Z_2^n��Ԫ�� Ĭ��Ϊ{-1,1}^n, type>=0��ʱ���� {0,1}^n
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
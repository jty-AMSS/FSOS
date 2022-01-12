function f=CZset(f,I,c)
%set \overwide(f)(T(i,:))=k(i);
A=f.c;
if size(I,2)==1
    A(I)=c;
    f.c=A;
else
    I=mod(I,2);
    I(I==0)=2;
    Index=zeros(size(I,1),1);
    for i=1:size(I,1)
        Index(i)=Z22Qindex(I(i,:));
    end
    %disp(Index)
     A(Index)=c;
     f.c=A;
end
end
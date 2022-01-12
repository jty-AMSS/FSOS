function v=CZget(f,I)
A=f.c;
if size(I,2)==1
    v=full(A(I));
else
    Index=zeros(size(I,1),1);
    for i=1:size(I,1)
        Index(i)=Z22Qindex(I(i,:));
    end
    %disp(Index)
     v=full(A(Index));
end
end

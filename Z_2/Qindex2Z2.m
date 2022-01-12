function I=Qindex2Z2(i,siz)
%Qindex to Z_2^n : N \mapsto n*1 vectors {1,2}^n 
% s.t.  T is a 2*2*2...*2 double array, T(ind2Z2n(siz,ndx))=T(ndx)

if length(siz)>1
    siz=length(siz);
end
T=dec2bin(i-1,siz);
I=(T(end:-1:1)-47);
end
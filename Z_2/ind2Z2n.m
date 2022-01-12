function I= ind2Z2n(siz,ndx)
%ind2Z2n : N \mapsto n*1 vectors {1,2}^n 
% s.t.  T is a 2*2*2...*2 double array, T(ind2Z2n(siz,ndx))=T(ndx)

if length(siz)==1
    siz=2+zeros(siz,1);
end
nout =length(siz);% max(nargout,1);
siz = double(siz);
lensiz = length(siz);

if lensiz < nout
    siz = [siz ones(1,nout-lensiz)];
elseif lensiz > nout
    siz = [siz(1:nout-1) prod(siz(nout:end))];
end

if nout > 2
    k = cumprod(siz);
    for i = nout:-1:3
        vi = rem(ndx-1, k(i-1)) + 1;
        vj = (ndx - vi)/k(i-1) + 1;
        var{i-2} = double(vj);
        ndx = vi;
    end
end

if nout >= 2
    vi = rem(ndx-1, siz(1)) + 1;
    v2 = double((ndx - vi)/siz(1) + 1);
    v1 = double(vi);
else 
    v1 = double(ndx);
end
I=[v1,v2,cell2mat(var)];
end
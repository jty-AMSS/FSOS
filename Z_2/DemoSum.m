%Demo sum(x[i],i=1,2,...,n)
n=10;
g=2-eye(n);g=[g;ones(1,n)*2];
c=ones(n,1);c=[c;n];
f=CZ_2n(g,c);
FF=CZifft(f);
min(FF(:))
t0=cputime;
F=ComputeSparseSOS(f);
Fsos=GenerateSOSForF(F);
norm(Fsos.c-f.c)
t=cputime-t0;
disp(['time:'])
disp(t);
x=sym('x',[1,n]);
for i=1:length(F)
    ff{i}=dispCZFunction(F{i});
end
for i=1:length(F)
    disp(vpa(ff{i}(x),5))
end
function dispCZValue(v)
if isstruct(v)
    v=CZifft(v);
end
N=size(v);
v=v(:);
for i=1:numel(v)
    t=Qindex2Z2(i,N);
    t(t==1)=-1;
    t(t==2)=1;
    disp(['[', num2str((t)) ']:' ,num2str(v(i))])
end
end
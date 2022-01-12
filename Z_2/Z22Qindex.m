function k=Z22Qindex(g)
% {1,2}^n \mapsto N
g=(g(1)<47)*47+g;
g=char(g);
g=g(end:-1:1);
k=bin2dec(g);
k=k+1;
end
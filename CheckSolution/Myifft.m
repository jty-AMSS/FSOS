function Fa=Myifft(c)
F=[c(end,:);c];F=F(1:end-1,:);F=[F(:,end),F];F=F(:,1:end-1);
Fb=fft2(F);
Fb=Fb(end:-1:1,end:-1:1);
Fa=Fb;
end
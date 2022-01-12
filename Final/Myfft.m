function W=Myfft(Fa)
Vsq=Fa;
VVsq=[Vsq(:,end),Vsq];
VVsq=VVsq(:,1:end-1);
VVsq=[VVsq(end,:);VVsq];
VVsq=VVsq(1:end-1,:);
W=ifft2(VVsq);
W(end:-1:1,end:-1:1)=W;
end
function [Q,Qindex]=L0_ifft(c,Qindex)
% It can be down only in Z_N \times Z_n
% % % [~,Qindex]=MyGenerateQFor2dim(c,Qindex);
% % % toc
GN=size(Qindex,1);
F=[c(end,:);c];F=F(1:end-1,:);F=[F(:,end),F];F=F(:,1:end-1);
Fb=fft2(F);
Fb=Fb(end:-1:1,end:-1:1);
Fa=Fb;
Vsq=(GN*Fa).^0.5;
VVsq=[Vsq(:,end),Vsq];
VVsq=VVsq(:,1:end-1);
VVsq=[VVsq(end,:);VVsq];
VVsq=VVsq(1:end-1,:);
W=ifft2(VVsq);
W(end:-1:1,end:-1:1)=W;
W=conj(W);
for i=1:GN
x=Qindex(i,1);y=Qindex(i,2);
u(i,1)=W(x,y);
end
Q=abs(u).^2;Q=Q(:);
%GN=size(Qindex,1);
% f=Myifft(c);
% f=real(f);
% f=f.^0.5;
% u=Myfft(f);
% for i=1:GN
%     Q(i)=u(Qindex(i,1),Qindex(i,2));
% end
% Q=Q(:);
end

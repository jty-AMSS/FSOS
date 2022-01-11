%TestAll
clear S;
clear Sc
clear GN
clear err;

S=0;
GN=numel(c{Index});
% for i=1:length(G{Index})
%     S=S+1/GN*abs(Myifft(G{Index}{i})).^2;
% end
% SS{Index}=S;
% Sc=Myfft(sym(S));
for i=1:length(G{Index})
    S=S+1/GN*ConjProd3(G{Index}{i},G{Index}{i});
end
SS{Index}=S;
Sc=S;
err=Sc-c{Index};
err=double(err);
err=abs(err(:));
err=max(err);
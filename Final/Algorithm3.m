function [Tp,chip]=Algorithm3(C,Qindex,T,chi,N)
% Input:Subsets  such that C{i}=Qindex(C_i,:);
% T:fourier support for example T={Qindex(1,:);Qindex(3,:);Qindex(4,:)}
% chi: chi[i]*C[i]  \subseteq T for all i \in length(C)
% chi=QIndex(Index of chi,:);
% output:local minimal Fourier support Tp
Tp=T;
chip=chi;
T=Tp;k=0;
for i=1:size(T,1)
    Tt=Tp;Tt(i-k,:)=[];
    [flag,chi]=CheckFourier(Tt,C,Qindex,N);
    if flag==1
        k=k+1;
        Tp=Tt;
        chip=chi;
    end
end
end

function [flag,chi]=CheckFourier(T,C,Qindex,N)
flag=1;
chi=zeros(size(C,1),size(Qindex,2));
for i=1:length(C)
    flagc=0;% not exist chi[i] :chi[i]*C[i]\subseteq  T
    for j=1:size(Qindex,1)% chi[j]*C[i]\subseteq  T
        Ct=GroupAdd(Qindex(j,:),C{i},N);
        Nt=size(unique([Ct;T],'rows'),1);
        Nc=size(T,1);
        if Nc==Nt % if Qindex[j]*C[i]\subseteq  T
            flagc=1;
            chi(i,:)=Qindex(j,:);
            break
        end
    end
    if flagc==0 % not exist chi[i] :chi[i]*C[i]\subseteq  T, T is not Fourier support
        flag=0;
        chi=[];
        break
    end
end
end


function z=GroupAdd(x,y,N)
%output x+y (in Z_N)
z=x+y;
N=ones(size(z,1),1)*N;
z=mod(z,N);
z(z==0)=N(z==0);
end
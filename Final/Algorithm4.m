function [Q,Gamma,Mvk,Nvk,chivk,T,Gammaalpha,Tp,chip,H]=Algorithm4(Q,Qindex,N)
% output Q is the 残差
%Prepartion Step
QI=find(diag(Q));
Q0=Q;
Q=Q(QI,QI);
Qindexp=Qindex(QI,:);
Gamma=(Q~=0);
Gamma=Gamma-diag(diag(Gamma));
Gammaalpha=Gamma;
Vnum=zeros(1,length(Q));Vnum=(Vnum==1);
T=[];
%step5
while(~all(Vnum)) % exists unnumbered
    v=find(Vnum==0,1);
    Vnum(v)=1;
    Nv{v}=(Gamma(v,:)==1);
    Nv{v}(v)=1;
    Gamma(Nv{v},Nv{v})=1;
    Gamma(v,:)=0;
    Gamma(:,v)=0;
    Nv{v}=Qindexp( Nv{v},:);
    %step9
    [chi{v},T]=ArgminStep9(T,Nv{v},Qindex,N);
    M{v}=GenerateMv(Q,v);
    if norm(M{v})==0
        M{v}=0;
    end
    Q=Q-M{v};
    Gamma=Gamma-diag(diag(Gamma));
    Gammaalpha=Gammaalpha+Gamma;
end
% step14;
Gammaalpha=Gammaalpha>0;
[nk,I]=InclusionMaximalElement(Nv);
for i=1:length(nk)
    Nvk{i}=Nv{nk(i)};
    chivk{i}=chi{nk(i)};
    Stemp=0;
    for j=1:length(I{i})
        Stemp=Stemp+M{I{i}(j)};
    end
    Mvk{i}=Stemp;
    if norm( Mvk{i})==0
        Mvk{i}=0;
    end
end
for i=1:length(Mvk)
    TempM=diag(Mvk{i});
    TempM=abs(TempM);
    MiIndex=find(TempM<1e-5);
    Mvk{i}(MiIndex,MiIndex)=0;
    TempM=diag(Mvk{i});
    TempM=abs(TempM);
    MiIndex=find(TempM>1e-5);
    TempM=Mvk{i}(MiIndex,MiIndex);
    TempM=1/2*(TempM+TempM');% M=M^*
    [L,D]=ldl(TempM);%L*D*L'=M{i};
    L=L*sqrt(max(D,0));
    H{i}=zeros(length(Mvk{i}));
    H{i}(MiIndex,MiIndex)=L';
    tempH = sparse(size(Qindex,1),size(Qindex,1));
    tempH(QI,QI)=H{i};
    H{i}=tempH;
end
[Tp,chip]=Algorithm3(Nvk,Qindex,T,chivk,N);
end

function z=GroupAdd(x,y,N)
%output x+y (in Z_N)
z=x+y;
N=ones(size(z,1),1)*N;
z=mod(z,N);
z(z==0)=N(z==0);
end

function [chiv,TuNv]=ArgminStep9(T,Nv,Qindex,N)
M=inf;chiv=[];
%for i=1:size(Qindex,1)
for i=size(Qindex,1):-1:1
    chi=Qindex(i,:);
    Tt=GroupAdd(chi,Nv,N);
    Tt=unique([T;Tt],'rows');
    if size(Tt,1)<M
        chiv=chi;
        M=size(Tt,1);
        TuNv=Tt;
    end
end
end


function Mv=GenerateMv(A,v)
n=length(A);
t=(A(v,:))~=0; % t=adj(v)
Mv=zeros(n,n);
if (abs(A(v,v))>1e-5)
    Mv(v,v)=A(v,v);
    Mv(v,t)=A(v,t);
    Mv(t,v)=A(t,v);
    Mv(t,t)=A(t,v)*A(t,v)'/A(v,v);
end
end

function [nk,I]=InclusionMaximalElement(Nv)
%输入：N{} cell array
%Output: index of maximal elemnts; eg nk=[1,5,7] then N{1},N{5},N{7} is the maximal elemnts
% I is the index that I{j} is the index subset of N{nk(j))  eg:
% nk=[1,5,7],I{2}=[2,4,5] then N{2},N{4},N{5} \subseteq  N{5}
n=length(Nv);
F=zeros(n);
for i=1:length(Nv)
    for j=1:length(Nv)
        C = intersect(Nv{i},Nv{j},'rows');
        F(i,j)=(size(C,1)==size(Nv{i},1));
        F(j,i)=(size(C,1)==size(Nv{j},1));
        if (F(i,j)==1&&F(j,i)==1)
            F(i,j)=0;
        end
    end
end
k=1;
FF=F;
for i=1:length(Nv)
    t=F(i,:);
    if sum(t)==0 %最大元
        nk(k)=i;
        t2=FF(:,i);
        I{k}=[find(t2)',i];
        FF(I{k},:)=0;
        k=k+1;
    end
end
end
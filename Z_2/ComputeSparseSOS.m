function [F,Qm]=ComputeSparseSOS(f)
%compute sparse SOS over Z_2^n
% input: f is a C[Z_2^n] structure
% output:F is a cell array such that f=sum(F{i}*F{i})
cvx_begin;cvx_end;cvx_clear;
disp('Start Compute')
n=f.n;
ti=cputime;
rootf=CZfft(CZifft(f).^0.5);
rc=full(rootf.c);
[T,y]=sort(abs(rc));
y=y(end:-1:1);
cvx_status='0';
MaxE=1;
MinE=1;
T=sort(T);
T=T(end:-1:1);
Tc=cumsum(T)./sum(T);
%Tcc=1-T./T(1);
tol=0.5;
GN=2^f.n;
%% 放大过程（找到大的MaxM使得 S_MaxM有解)
 cvx_precision high
while(~strcmp(cvx_status,'Solved'))% if is not satisified
    k=find(Tc>tol,1);
    if(k<=MinE)
        k=MinE+1;
    end
    %disp(k)
    %disp(1-tol)
    GNp=k;
    %增加波及到的B
    g=[];
    for i=1:k
        for j=1:k
            g=[g;GroupAdd(Qindex2Z2(y(i),n),-Qindex2Z2(y(j),n))];
        end
    end
    g1=find(f.c);
    g2=[];
    for i=1:length(g1)
        g2=[g2;Qindex2Z2(g1(i),n)];
    end
    g=[g;g2];g=unique(g,'rows');
    clear Bp;
    for i=1:size(g,1)
        Bp{i}=GenerateBchiPrime2(n,g(i,:),y(1:k));
        Bp{i}=Bp{i}(y(1:k),y(1:k));
    end
    cvx_solver Mosek
    cvx_begin sdp quiet
    variable Q(GNp,GNp) hermitian
    minimize(1)
    for i=1:size(g,1)
        Bp{i}(:)'*Q(:)==CZget(f,g(i,:));
    end
    Q >= 0
    cvx_end
    if(strcmp(cvx_status,'Solved'))
        MaxE=k;
        Qm=sparse(GN,GN);
        Qm(y(1:k),y(1:k))=Q;
    else
        tol=1-tol;
        tol=tol*0.5;
        tol=1-tol;
        MinE=k;
    end
end
%% Searching
while(MaxE-MinE>1) 
    k=round(0.5*(MaxE+MinE));
    GNp=k;
    g=[];
    for i=1:k
        for j=1:k
            g=[g;GroupAdd(Qindex2Z2(y(i),n),-Qindex2Z2(y(j),n))];
        end
    end
    g1=find(f.c);
    g2=[];
    for i=1:length(g1)
        g2=[g2;Qindex2Z2(g1(i),n)];
    end
    g=[g;g2];g=unique(g,'rows');
    clear Bp;
    for i=1:size(g,1)
        Bp{i}=GenerateBchiPrime2(n,g(i,:),y(1:k));
        Bp{i}=Bp{i}(y(1:k),y(1:k));
    end
    cvx_begin sdp quiet
    variable Q(GNp,GNp) hermitian
    minimize(1)
    %LabelOfChi0InY=find(y==GN);
    %minimize(-Q(LabelOfChi0InY,LabelOfChi0InY))
    for i=1:length(g)
        Bp{i}(:)'*Q(:)==CZget(f,g(i,:));
    end
    Q >= 0
    cvx_end
    if(strcmp(cvx_status,'Solved'))
        MaxE=k;
        Qm=sparse(GN,GN);
        Qm(y(1:k),y(1:k))=Q;
    else
        MinE=k;
    end
end
%toc
disp(['sparsity:' num2str(length(find(diag(Qm))))])

y=find(diag(Qm));
[L,D]=ldl(full(Qm(y,y)));%L*D*L'-full(F)
L=2^(-f.n/2)*sqrt(D)*L';
gIndex=zeros(length(y),f.n);
for i=1:length(y)
    gIndex(i,:)=Qindex2Z2(y(i),f.n);
end
F=cell(length(y),1);
for i=1:size(L,1)
    F{i}=CZ_2n(gIndex,L(i,:).');
end
end

function z=GroupAdd(x,y)
z=x+y;
z=mod(z,2);
z(z==0)=2;
end


function B=GenerateBchiPrime2(N,g,I)
% Bp{i}=GenerateBchiPrime2(n,g(i,:),y(1:k))
% Generate The B_chi for
% chi=g \in Z_2^N
%I= Index in Qindex such that B_chi only not vanish at  label I, e.g.
%I=[1;4;5;6], mens B is the submatrix of B_chi(I,I)
Qindex=zeros(length(I),N);
for i=1:length(I)
    Qindex(i,:)=Qindex2Z2(I(i),N);
end
J=GroupAdd(Qindex,g);
J=intersect(J,Qindex,'rows');
Jm=GroupAdd(J,-g);%B(j,jm)=1/GN
ind1=zeros(size(J,1),1);
ind2=ind1;
for i=1:size(J,1)
    %     ind1(i)=find(all(Qindex'==Jm(i,:)'));
    %     ind2(i)=find(all(Qindex'==J(i,:)'));
    ind1(i)=Z22Qindex(Jm(i,:));
    ind2(i)=Z22Qindex(J(i,:));
end
GN=2^N;
B=sparse(ind1,ind2,1/GN,GN,GN);
end
function [Qm,MaxE,T]=AfterL0SDP3(Qg,c,Qindex)
% It can be down only in Z_N \times Z_n
% 计算精化结果
% if method=1 we use 二分 =0 we add again and again
% if method<1 then we check S_k,S_{k+1}... until S_{k+m} is not empty ,here
% k=#{Q(x,x):Q(x,x)>method}
GN=length(c(:));
T=abs(Qg);
[~,y]=sort(T);
y=y(end:-1:1);
cvx_status='0';
MaxE=1;
MinE=1;
T=sort(T);
T=T(end:-1:1);
Tc=cumsum(T)./sum(T);
Tcc=1-T./T(1);
tol=0.9;
%% Search phase
while(~strcmp(cvx_status,'Solved'))% if is not satisified
    k=find(Tc>tol,1);
    while(k==MinE)
        tol=1-tol;
        tol=tol*0.5;
        tol=1-tol;
        k=find(Tcc>tol,1);
    end
    GNp=k;
    g=[];
    for i=1:k
        for j=1:k
            g=[g;GroupAdd(Qindex(y(i),:),-Qindex(y(j),:),size(c))];
        end
    end
    [g1,g2]=find(c);g1=g1(:);g2=g2(:);
    g=[g;[g1,g2]];g=unique(g,'rows');
    clear Bp;
    for i=1:size(g,1)
        Bp{i}=GenerateBchiPrime2(size(c),g(i,:),Qindex,y(1:k));
        Bp{i}=Bp{i}(y(1:k),y(1:k));
    end
    cvx_begin sdp quiet
    variable Q(GNp,GNp) hermitian
    minimize(1)
    for i=1:length(g)
        Bp{i}(:)'*Q(:)==c(g(i,1),g(i,2));
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
%% binary search  phase
while(MaxE-MinE>1) 
    k=round(0.5*(MaxE+MinE));
    GNp=k;
    g=[];
    for i=1:k
        for j=1:k
            g=[g;GroupAdd(Qindex(y(i),:),-Qindex(y(j),:),size(c))];
        end
    end
    [g1,g2]=find(c);g1=g1(:);g2=g2(:);
    g=[g;[g1,g2]];g=unique(g,'rows');
    clear Bp;
    for i=1:size(g,1)
        Bp{i}=GenerateBchiPrime2(size(c),g(i,:),Qindex,y(1:k));
        Bp{i}=Bp{i}(y(1:k),y(1:k));
    end
    %cvx_solver SDPT3
    cvx_begin sdp quiet
    variable Q(GNp,GNp) hermitian
    %minimize(1)
    LabelOfChi0InY=find(y==GN);
    minimize(-Q(LabelOfChi0InY,LabelOfChi0InY))
    for i=1:length(g)
        Bp{i}(:)'*Q(:)==c(g(i,1),g(i,2));
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
disp(['sparsity:' num2str(length(find(diag(Qm))))])
end

%% Tool Function


function B=GenerateBchiPrime2(N,g,Qindex,I)
% Generate The B_chi for
% chi=Qindex(g,:) or g is just the chi (size(g)=[1,2])
% N=size(c)
%I= Index in Qindex such that B_chi only not vanish at  label I, e.g.
%I=[1;4;5;6], mens B is the submatrix of B_chi(I,I)
if length(g)==1
    g=Qindex(g,:);
end
J=GroupAdd(Qindex(I,:),g,N);
J= intersect(J,Qindex(I,:),'rows');
Jm=GroupAdd(J,-g,N);%B(j,jm)=1/GN
ind1=zeros(size(J,1),1);
ind2=ind1;
for i=1:size(J,1)
    ind1(i)=find(all(Qindex'==Jm(i,:)'));
    ind2(i)=find(all(Qindex'==J(i,:)'));
end
GN=prod(N);
B=sparse(ind1,ind2,1/GN,GN,GN);
end

function z=GroupAdd(x,y,N)
%output x+y (in Z_N)
z=x+y;
z=mod(z,N);
z(z(:,1)==0,1)=N(1);
z(z(:,2)==0,2)=N(2);
end


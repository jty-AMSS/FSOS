function [G,AllG]=GenerateSparseSoS(chip,H,Qindex,N)
%chip{i}=Qindex(j,:) s.t. \chi{j}*Nv{j} \subseteq T
% H{i}'*H{i}=M{i};
%i=1,2,...,|G|.
k=1;
for i=1:length(H)
    %    if norm(H{i})~=0 %H{i}!=0
    if full(max(abs(H{i}(:))))~=0
        for j=1:length(H{i})
            if norm(H{i}(j,:))~=0 %H{i}(j,:)!=0
                 G{k}=zeros(N);
                %G{k}=sparse(N(1),N(2));
                temp=H{i}(j,:);
                Index=find(temp);
                for t=Index
                    x=Qindex(t,:);
                    z=GroupAdd(x,chip{i},N);
                    G{k}(z(1),z(2))=H{i}(j,t);
                end
                k=k+1;
            end
        end
    end
end
AllG=zeros(N);
for i=1:length(G)
    AllG=AllG+(abs(G{i})~=0);
end
% VarChi=GenerateChi(N);
% g=sym(0);
% for i=1:length(G)
%     tempg=sym(G{i}).*VarChi;
%     tempg=sum(tempg(:));
%     g=g+1/prod(N)*abs(tempg)^2;
% end
end


function z=GroupAdd(x,y,N)
%output x+y (in Z_N)
z=x+y;
N=ones(size(z,1),1)*N;
z=mod(z,N);
z(z==0)=N(z==0);
end
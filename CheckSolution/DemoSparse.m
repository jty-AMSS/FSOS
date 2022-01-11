%DemoSparse:
%Second experiment: bounded support and bounded minimum
%Please call this program in "MeanDeom.m"
fSparsity=12; %  Bound of sparsity of f, |supp(f)|=2*fSparsity+1
cvx_precision  DEFAULT
for Index=1:1
    cvx_begin;cvx_end;cvx_clear;cvx_solver Mosek
    n=(50*k).^2;m=1; %Group :Z_n \times Z_m
    cp=zeros(n,m);N=[n,m];
    for i=1:fSparsity
        t=2*(rand()-1/2)+1i*2*(rand()-1/2);
        x=randi(n);y=randi(m);
        cp(x,y)=cp(x,y)+t;
        z=GroupAdd(-[x,y],0,N);
        cp(z(1),z(2))=cp(z(1),z(2))+t';
    end
    cp(end)=0;
    c{Index}=cp;
    V=Myifft(cp);
    cp(end)=-min(real(V(:)))+rand();
    c{Index}=cp/max(abs(cp(:)));
    cd ..;
    cd Final
    DemoPure
    cd ..
    cd CheckSolution
    disp(length(unique([x,y],'rows')))
    %err=err*max(abs(cp(:)));
    str{Index}=(['sparsity of f: ',num2str(length(find(c{Index}))), '.Local minimal Fourier suppiorts: ' num2str(size(Tp{Index},1)) , ' with coefficient error: ' num2str(err),[' Time: ',num2str(Tim(Index))]]);
    str2{Index}=['$ \mathbb{Z}_{', num2str(n), '}\times \mathbb{Z}_{', num2str(m), '}$ & ',num2str(length(find(c{Index}))) ,'&', num2str(size(Tp{Index},1))  ,' &Mosek& ',num2str(Tim(Index)),'&',num2str(err) ,'\\'];
end


function z=GroupAdd(x,y,N)
%output x+y (in Z_N)
z=x+y;
z=mod(z,N);
z(z(:,1)==0,1)=N(1);
z(z(:,2)==0,2)=N(2);
end
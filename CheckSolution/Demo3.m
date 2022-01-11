% Third experiment: bounded support
k=7;s=15;
cvx_precision DEFAULT
%cvx_precision high
for Index=1:5
    cvx_begin;cvx_end;cvx_clear;cvx_solver Mosek
    m=(Index*50);n=(50*Index);
    nzx=randi(n);nzy=randi(m);
    N=[n,m];
    cp=zeros(n,m);
    G0{Index}=zeros(n,m);
    G1{Index}=zeros(n,m);
    for i=1:k
        x=randi(n);y=randi(m);
        if ~(x==n&&y==m)
            G0{Index}(x,y)=1;
        end
    end
    for i=1:s
        z=G0{Index}.*(20*(rand(n,m)-1/2))+G0{Index}.*1i.*(20*(rand(n,m)-1/2));
        Z{i}=z;
        cp=cp+ConjProd3(z,z);
    end
    c{Index}=cp;
    V=Myifft(cp);
    %disp mean value of $f$:
    min(real(V(:)))
    c{Index}=cp/max(abs(cp(:)));
    cd ..;
    cd Final
    DemoPure
    cd ..
    cd CheckSolution
    %disp(length(unique([x,y],'rows')))
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
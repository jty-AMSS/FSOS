%DemoDegree:
%First experiment:bounded degree and bounded minimum
%Please call this program in "MeanDeom.m"
deg=24;%Degree bound of f
for Index=1:1
    cvx_begin;cvx_end;cvx_clear
    cvx_solver Mosek
    cvx_precision DEFAULT
    n=(50*k)^2;m=1;
    %n=(50*k);m=(50*k);
    N=[n,m]; %Group :Z_n \times Z_m
    cp=zeros(n,m);
    dxm=min(deg,floor(n/2));
    dym=min(deg,floor(m/2));
    for i=-dxm:dxm
        for j=-dym:dym
            x=mod(i,n);
            y=mod(j,m);
            z=GroupAdd([x,y],[0,0],N);
            x=z(1);y=z(2);
            %t=20*(rand()+1i*rand())-10-10*1i;
            t=2*(rand()-1/2)+1i*2*(rand()-1/2);
            cp(x,y)=cp(x,y)+t;
            z=GroupAdd(-[x,y],[0,0],N);
            x=z(1);y=z(2);
            cp(x,y)=+cp(x,y)+t';
        end
    end
    cp(end)=0;
    V=Myifft(cp);
    cp(end)=-min(real(V(:)))+rand();
    c{Index}=cp/max(abs(cp(:)));
    clear g;
    cd ..;
    cd Final
    DemoPure
    cd ..
    cd CheckSolution
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

%DemoPure
disp('Start Compute')
Qindex=GenerateQindex(c{Index});
ti=cputime;
[Qg{Index},Qindex]=L0_ifft(c{Index},Qindex);
[Qp{Index},MaxM{Index},ST{Index}]=AfterL0SDP3(Qg{Index},c{Index},Qindex);% with Ã¿´Î+1
% allow Numberic error:
QpI=diag(Qp{Index});QpI=abs(QpI)>0;
Qp2=sparse(size(Qp{Index},1),size(Qp{Index},2));
Qp2(QpI,QpI)=Qp{Index}(QpI,QpI);
[Q{Index},Gamma{Index},Mvk{Index},Nvk{Index},chivk{Index},T{Index},Gammaalpha{Index},Tp{Index},chip{Index},H{Index}]=Algorithm4(Qp2,Qindex,N);
Tim(Index)=cputime-ti;
disp('Group:')
disp([n,m])
disp('Computation ended , now Check Solution')
cd ..
cd CheckSolution
[Gp,AllGp]=GenerateSparseSoS(chip{Index},H{Index},Qindex,N);
G{Index}=Gp;
AllG{Index}=AllGp;
TestAll
disp(['sparsity of f: ',num2str(length(find(c{Index}))), '.Local minimal Fourier suppiorts: ' num2str(size(Tp{Index},1)) , ' with coefficient error: ' num2str(err)])
disp(['Time: ',num2str(Tim(Index))])
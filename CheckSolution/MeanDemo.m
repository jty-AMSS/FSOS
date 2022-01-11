%Mean Demo: Test many times and output the mean value of experiments
for k=1:5 %k discripts the size of group
for k2=1:10 % k2 is the number of replication of experiments
t1=cputime;
%DemoDegree; %First experiment:bounded degree and bounded minimum
DemoSparse; %Second experiment: bounded support and bounded minimum
t2=cputime;
Tim2(k,k2)=t2-t1;%Time costs
sparsity(k,k2)=(length(find(c{Index})));%sparsity
LMFS(k,k2)=(size(Tp{Index},1));%Local minimal Fourier suppiorts
Err(k,k2)=err; % error
end
end

%% Z_N *Z_N
% for Index=1:k
%     str2{Index}=['$ \mathbb{Z}_{', num2str(Index*20), '}\times \mathbb{Z}_{', num2str(Index*20), '}$ & ',num2str(mean(sparsity(Index,:))) ,'&', num2str(mean(LMFS(Index,:)))  ,' &Mosek& ', num2str(mean(Tim2(Index,:))) ,'&',num2str(mean(Err(Index,:))) ,'\\'];
% end
%% Z_N
for Index=1:k
    str2{Index}=['$ \mathbb{Z}_{', num2str(Index^2*2500), '}$ & ',num2str(mean(sparsity(Index,:))) ,'&', num2str(mean(LMFS(Index,:)))  ,' &Mosek& ', num2str(mean(Tim2(Index,:))) ,'&',num2str(mean(Err(Index,:))) ,'\\'];
end
str3=str2;
save str2
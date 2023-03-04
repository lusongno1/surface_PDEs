%convergence rate plot for eigenvalues
clc
clear
close all
%%导入数据
k = 2;
%data = [];
j = 1;
%N_sets = [8 16 32 64 128];
N_sets = [8 16 32 64 128];
Nk = size(N_sets,2);
for N = N_sets
    %path = ['NISO_k' num2str(k) 'N' num2str(N)];
    path = ['ngsxfem_stb_k' num2str(k) 'N' num2str(N)];
    Eigpath = [path '/Eig_Clu.mat'];
    data(j) = load(Eigpath);
    j = j+1;
end

for s=1:10
    for i=1:length(data)
        Eig_Clu = data(i);
        Eig_Clu = Eig_Clu.Eig_Clu;
        realEig = Eig_Clu{s,1};
        calEig = Eig_Clu{s,2};        
        %Err(s,i) = norm( abs(realEig-calEig),1);
        Err(s,i) = sum(abs(realEig-calEig))/length(calEig);
        %Err(s,i) = abs(Err(s,i)/realEig);
    end
end
legendname = [];

one_dx = N_sets;
dx = 1./N_sets;
for i = 1:10
    %semilogy(dx,Err(i,:),'-o')
    %semilogy(1./dx,Err(i,:),'-o')  
    pic(i) = loglog(one_dx,Err(i,:),'-o')
    %plot(log10(1./dx),log10(Err(i,:)),'-o')
    %plot(log2(1./dx),log2(Err(i,:)),'-o')
    %plot(log10(Err(i,:)./realEig),'-o')    
    hold on;
    %legendname(end+1) =  ['Lam' num2str(i)];
end

for i = 1:10
log2(Err(i,1:Nk-1)./Err(i,2:Nk))
end

lgd = legend;

fig_beauty;
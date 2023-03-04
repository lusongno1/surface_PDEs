%a small test about error calculation
clc
clear
close all
format long
for i = 8%i=[8 16 32 64 128]
    path = ['NISO_k1N' num2str(i)] 
    Eigpath = [path '/Eig_Clu.mat'];
    Eig = load(Eigpath);
    Eig_Clu = Eig.Eig_Clu;
    
    for i=1:length(Eig_Clu)
        len(i) = length(Eig_Clu{i,2});
        realEig = Eig_Clu{i,1};
        calEig = Eig_Clu{i,2};
        Err(i) = sum(abs(realEig-calEig))/length(calEig);
        Err(i) = abs(Err(i)/realEig);
    end
    len
    Err
end
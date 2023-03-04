%make a zero padding for eig_cal and Eig_Clu
clc
clear
close all
%%导入数据
path = 'NISO_k2N128';
Apath = [path '/A.txt'];
Bpath = [path '/M.txt'];
load(Apath)
load(Bpath)
A = sparse(A(:,1),A(:,2),A(:,3));
B = sparse(M(:,1),M(:,2),M(:,3));
eig_cal_fix = eigs(A,B,1,0);
load([path '/eig_cal.txt']);
if(eig_cal(1)>0.001)
    eig_cal = [eig_cal_fix;eig_cal];
end
dlmwrite([path '/eig_cal.txt'],eig_cal,'delimiter','\n','precision',15)
%%
k = 1:10;
eig_true = k.*(k-1);
Eig_Clu = eig_cluster2(eig_true,eig_cal);
Eig_Clu = Eig_Clu(1:10,:);
save([path '/Eig_Clu.mat'],'Eig_Clu')

fileID = fopen([path '/Eig_Clu.csv'],'w');
fprintf(fileID,'%s,%s\n',"Real Eigenvalue","Numerical Eigenvalues");
for i=1:length(Eig_Clu)
    realEig = Eig_Clu{i,1};
    calEigs = Eig_Clu{i,2};
    fprintf(fileID,'%d,',realEig);
    for j=1:length(calEigs)
        fprintf(fileID,'%.15f ',calEigs(j)); 
    end
    fprintf(fileID,'\n'); 
end
fclose(fileID);
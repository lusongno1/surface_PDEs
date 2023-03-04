%calculate eigenvalues by method of singgep from rcp paper
clc
clear
close all
%%
path = 'NISO_k2N32';
Apath = [path '/A.txt'];
Bpath = [path '/M.txt'];
load(Apath)
load(Bpath)
if exist([path '/log.txt'],'file')==2
    delete([path '/log.txt']);
end
diary([path '/log.txt']);
diary on;
A = sparse(A(:,1),A(:,2),A(:,3));
B = sparse(M(:,1),M(:,2),M(:,3));
size = size(A,1);
%%
A = full(A);
B = full(B);
[lambda,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,1);
%[lambda,Z,d,X,Y,U,V,DA,DB] = singgep_modified(A,B,1);
lambda = sort(lambda);
eig_cal = lambda;
dlmwrite([path '/eig_cal.txt'],eig_cal,'delimiter','\n','precision',15)
%%
k = 1:10;
eig_true = k.*(k-1);
%eig_true = [0 eig_true];
Eig_Clu_origin = eig_cluster2(eig_true,eig_cal);
Eig_Clu = Eig_Clu_origin(1:10,:);
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


diary off;

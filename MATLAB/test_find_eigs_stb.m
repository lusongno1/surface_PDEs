%find eigenvalues with A and B by MATLAB function eig and collocate eigenvalues as result of eig_cluster
%add stabilization term
clc
clear
close all

%% 
%path = 'ISO_k1N32_Stb';
%path = 'NISO_k2N32';
path = 'ngsxfem_stb_k2N128';
Apath = [path '/A.txt'];
Bpath = [path '/M.txt'];
Spath = [path '/S.txt'];
load(Apath)
load(Bpath)
load(Spath)
if exist([path '/log.txt'],'file')==2
    delete([path '/log.txt']);
end
diary([path '/log.txt']);
diary on;
%% 
A = sparse(A(:,1),A(:,2),A(:,3));
B = sparse(M(:,1),M(:,2),M(:,3));
S = sparse(S(:,1),S(:,2),S(:,3));
%updata A B and S
sz = size(S,1);
sizeA = size(A,1);
A = A+S;

%% cal calculation of eigenvalues
%eig_cal = eigs(A,B,100,'smallestabs');
%eig_cal = real(eig_cal);
eig_cal = eigs(A,B,150,-0.1);

%A_full = full(A);
%B_full = full(B);
%% 
dlmwrite([path '/eig_cal.txt'],eig_cal,'delimiter','\n','precision',15)
%% 
k = 1:11;
eig_true = k.*(k+1);
eig_true = [0 eig_true];
Eig_Clu_origin = eig_cluster(eig_true,eig_cal);
Eig_Clu = Eig_Clu_origin(1:11,:);
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
%% 
% scatter(eigenvalues,zeros(1,length(eigenvalues)))
% hold on
% scatter(eigenvalues_cal,zeros(1,length(eigenvalues_cal)))



%%

%plot(eigenvalues_cal)
%hold on
%plot(eigenvalues)

%%
for i=1:length(Eig_Clu)
    scatter(i,Eig_Clu{i,1},'red','filled');
    hold on;
    calEigs = Eig_Clu{i,2};
    scatter(repmat(i,1,length(calEigs)),calEigs,'green');
end

legend('Real eigenvalue','Numerical eigenvalues');

%%

for s=1:10
    realEig = Eig_Clu{s,1};
    calEig = Eig_Clu{s,2};
    %Err(s,i) = norm( abs(realEig-calEig),1);
    Err(s) = sum(abs(realEig-calEig))/length(calEig);
end
format long;
%Err
%vpa(Err,10)



diary off;

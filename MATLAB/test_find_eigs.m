%find eigenvalues with A and B by MATLAB function eig，collocate eigenvalues as result of eig_cluster
clc
clear
close all

%% 导入数据
path = 'NISO_k1N32_Stb';
%path = 'NISO_k2N32';
Apath = [path '/A.txt'];
Bpath = [path '/M.txt'];
load(Apath)
load(Bpath)
if exist([path '/log.txt'],'file')==2
    delete([path '/log.txt']);
end
diary([path '/log.txt']);
diary on;
%% 生成 sparse 矩阵
A = sparse(A(:,1),A(:,2),A(:,3));
B = sparse(M(:,1),M(:,2),M(:,3));
size = size(A,1);

% A=sparse(A(:,1),A(:,2),A(:,3));
% B=sparse(M(:,1),M(:,2),M(:,3));
% size=size(A,1);


%eigA = eig(A);
%eigB = eig(B);

B_epsilon = sparse(1:size,1:size,1e-15);
B = B+B_epsilon;
%A = A+B_epsilon;
%B = B+diag(ones(size(B,1),1)*1e-15);%�? B 的对角元加个扰动

A_full = full(A);
B_full = full(B);

%size = size(A,1)%矩阵的规�?
%% 观察矩阵的数据分布，判断其是否对称正定，是否满秩
%spy(A)
% A - A.'
eigA = eig(A);
%rank(A_full)
rankA = sum(eigA>0)
isDefiniteA = ~sum(~(eigA >= 0))
%isDefiniteA = ~sum(~(eigA >= -1e-10))
%
% %spy(B)
% B - B.'
eigB = eig(B);
rankB = sum(eigB>0)
%condB = cond(B)
%cond(A)
%rank(B_full)
%isDefiniteB = ~sum(~(eigB >= -1e-10))
isDefiniteB = ~sum(~(eigB >= 0))
if ~isDefiniteB
    error('B is not definite!!!');
end

%% 求广义特征�?�问�?
%eigs(A,B,60,'smallestabs')
eig_cal = eig(A_full,B_full).';
%max(imag(eig_cal))
%[eig_vec,eig_cal] = eig(A_full,B_full);
dlmwrite([path '/eig_cal.txt'],eig_cal,'delimiter','\n','precision',15)
%dlmwrite([path '/eig_vec.txt'],eig_vec.','delimiter',',')%,'precision',5)

%eigenvalues_cal_unique = unique(eigenvalues_cal)
%eigenvalues_cal = eigenvalues_cal(1:100)
%imag(eig_cal)
%real(eig_cal)
%% 实际特征�?
k = 1:11;
eig_true = k.*(k+1);
eig_true = [0 eig_true];
Eig_Clu_origin = eig_cluster(eig_true,eig_cal);
Eig_Clu = Eig_Clu_origin(2:11,:);
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
%% 绘图
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

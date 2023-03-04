%calculate eigenvalues by MATLAB function eigs when data dimension is too large
clc
clear
close all
%%导入数据
path = 'ISO_k1N128';
Apath = [path '/A.txt'];
Bpath = [path '/M.txt'];
load(Apath)
load(Bpath)
if exist([path '/log.txt'],'file')==2
    delete([path '/log.txt']);
end
diary([path '/log.txt']);
diary on;
%% 生成矩阵
A = sparse(A(:,1),A(:,2),A(:,3));
B = sparse(M(:,1),M(:,2),M(:,3));
size = size(A,1);
%B_epsilon = sparse(1:size,1:size,1e-5);
%B = B+B_epsilon;
%B = B+diag(ones(size(B,1),1)*1e-15);
%A_full = full(A);
%B_full = full(B);

% %% 观察矩阵的数据分布，判断其是否对称正定，是否满秩
% %spy(A)
% % A - A.'
%eigA = eigs(A,6,0);
% 
% 
% 
% rank(A_full)
% 
%isDefiniteA = ~sum(~(eigA >= 0))
% % 
% % %spy(B)
% % B - B.'
%eigB = eigs(B,6,0);
% rank(B_full)
%isDefiniteB = ~sum(~(eigB >= 0))
%if ~isDefiniteB
%    error('B is not definite!!!');
%end
%% 实际特征值
k = 1:10;
eig_true = k.*(k-1);
%tmp = eigs(A,B,100,2)
%% 求广义特征值问题
%eig_cal = [];
% for i=1:length(eig_true)
%     i
%     sigma = eig_true(i);
%     eig_cal = [eig_cal;eigs(A,B,2*i+1,sigma)];
% end
eig_cal = eigs(A,B,150,0);
%eig_cal = eig_cal(eig_cal>-1);


dlmwrite([path '/eig_cal.txt'],eig_cal,'delimiter','\n','precision',15)

%%
%Eig_Clu = eig_cluster2(eig_true,eig_cal);
%imag(eig_cal)
%eig_cal = eig(A_full,B_full).';
%eigenvalues_cal_unique = unique(eigenvalues_cal)
%eigenvalues_cal = eigenvalues_cal(1:100)


%%
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
%% 绘图
% scatter(eigenvalues,zeros(1,length(eigenvalues)))
% hold on
% scatter(eigenvalues_cal,zeros(1,length(eigenvalues_cal)))



%%

%plot(eigenvalues_cal)
%hold on
%plot(eigenvalues)

%%
% for i=1:length(Eig_Clu)
%     scatter(i,Eig_Clu{i,1},'red','filled');
%     hold on;
%     calEigs = Eig_Clu{i,2};
%     scatter(repmat(i,1,length(calEigs)),calEigs,'green');
% end
% 
% legend('Real eigenvalue','Numerical eigenvalues')

diary off;


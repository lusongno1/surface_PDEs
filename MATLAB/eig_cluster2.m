function R = eig_cluster2(eig_true,eig_cal)% not efficient,but easy to implement
%eig_cal(eig_cal < 1e-6) = [];%过滤掉 0 特征值
for i=1:length(eig_true)
    begin = (i-1)^2+1;
    finish = i^2;
    R{i,1} =  eig_true(i);
    R{i,2} = eig_cal(begin:finish);
end

% for i=1:length(eig_cal)
%     [~,ind] = min(abs(eig_cal(i) - eig_true));
%     R{ind,2} = [R{ind,2} eig_cal(i)];
% end
    
end
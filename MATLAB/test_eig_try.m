% clc
% clear
% close all
% A  = [2 -1;-1 2];
% M = [2/3 1/6;1/6 2/3];
% A_= [-2 -1 -1;-1 2 2;-1 2 2];
% M_ = [2/3 1/6 1/6;1/6 2/3 2/3;1/6 2/3 2/3];
% M_ = M_+diag(ones(size(M_,1),1)*1e-10)
% 
% eig(A,M)
% eig(A_,M_,'qz')
% eig(A_,M_,'chol')
% eig(A_,M_,'vector')
% eigs(A_,M_,3)
% 
% eig(A)
% eig(M)
% [V,D] = eig(A_,M_)
% 
% [V,D] = eig(A_)
% [V,D] = eig(M_)
% 
% eig(A_)
% eig(M_)

clc
clear
A_in = diag([1 2 -1e-15])
%A_in = diag([1 2 0])
B_in = diag([1 1 1e-15])
A = A_in;
B = B_in;
%A = A+diag(ones(size(B,1),1)*1e-10)
%B = B+diag(ones(size(B,1),1)*1e-15)
[V,D]  = eig(A,B)
eig_cal = diag(D);
err = A_in*V - B_in*V*D;
max(  (abs(err)))





%a small test for calculating eigenvalues by data of ISO
clc
clear
close all

%% 
path = 'ISO_k2N8';
Apath = [path '/A.txt'];
Bpath = [path '/M.txt'];
load(Apath)
load(Bpath)
diary on;
%% 
A_in = sparse(A(:,1),A(:,2),A(:,3));
B_in = sparse(M(:,1),M(:,2),M(:,3));
A = A_in;
B = B_in;
N = size(A,1);

B_epsilon = sparse(1:N,1:N,1e-15);
A = A+B_epsilon;
B = B+B_epsilon;
%P = diag(diag(inv(B)));
%A = P*A;
%B = P*B;




A_full = full(A);
B_full = full(B);
%eigA = eig(A);

%A_full = balance(A_full);
%B_full = balance(B_full);

%% 
[V,D]  = eig(A_full,B_full,'chol');
%[V,D]  = eig(A_full,B_full,'qz');
eig_cal = diag(D);
%eig_cal =  real(eig_cal);
%err = A_full*V - B_full*V*D;
err = A_in*V - B_in*V*D;
max(max(abs(err)))

dlmwrite([path '/eig_cal.txt'],eig_cal,'delimiter','\n','precision',15)
%%
diary off;

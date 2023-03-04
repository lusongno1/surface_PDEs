%a small test about eigenvalues distributions
clc
clear
a = 1/3;
b = 1/6;
N = 5;
A=diag(repmat([2], 1, N))+diag(repmat([-1], 1, N-1), 1)+diag(repmat([-1], 1, N-1), -1);

A(1,end) = b;
A(end,1) =  b;

A

eig(A)



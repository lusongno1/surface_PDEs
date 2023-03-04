%an useless test about left eigenvalues and right eigenvalues
clc
clear
close all
N = 30
B = magic(N)
A = randi(N,N)
A = A+A.';
B = B+B';

D = eigs(A,5,'smallestabs')
[X,D,Y] = eig(A)
eigs(A,3,'smallestabs')

[X,D,Y] = eig(A,B)
[V1,D1] = eigs(A,B,N,'smallestabs')
[V2,D2] = eigs(A.',B.',N,'smallestabs')

diag(D1)-diag(D2)
%find eigenvalues by method singgep, just test
clc
clear
close all

%% 
path = 'NISO_k2N4';
Apath = [path '/A.txt'];
Bpath = [path '/M.txt'];
load(Apath)
load(Bpath)
%% 
A = sparse(A(:,1),A(:,2),A(:,3));
B = sparse(M(:,1),M(:,2),M(:,3));
A = full(A);
B = full(B);
%eigB = eig(B)
% A = [ -1    -1    -1    -1    -1    -1    -1
%        1     0     0     0     0     0     0
%        1     2     1     1     1     1     1
%        1     2     3     3     3     3     3
%        1     2     3     2     2     2     2
%        1     2     3     4     3     3     3
%        1     2     3     4     5     5     4];
%    
% B = [ -2    -2    -2    -2    -2    -2    -2
%        2    -1    -1    -1    -1    -1    -1
%        2     5     5     5     5     5     5
%        2     5     5     4     4     4     4
%        2     5     5     6     5     5     5
%        2     5     5     6     7     7     7
%        2     5     5     6     7     6     6];

[lambda,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,1);
%[lambda,Z,d,X,Y,U,V,DA,DB] = singgeps(A,B,2,1);
%l2 = eig(A,B)


lambda = sort(lambda)
%l2(2:end) - lambda
%max(imag(lambda))

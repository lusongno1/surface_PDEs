clc
clear
close all
A = diag([1,2,3,0,0,0])
B = diag([3,4,0,0,0,0])
[lambda,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,1);
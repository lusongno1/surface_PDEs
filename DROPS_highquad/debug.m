 clc
 clear
 %A1 = 0.24085*12+0.616645*12
 %A2 = 0.261799*12+0.785399*12
  %A1 = 0.252586*12+0.41363*12
 %A2 = 0.261799*12+0.0850216*12
 %%A = 4*pi;
 %%N = 16;
 %%for N=[2,4,8,16,32,64]
 %A1 = load(["data/N" num2str(N) "Origin.txt"]);
 %A2 = load(["data/N" num2str(N) "High.txt"]);
 %%AC = load(["new1/N" num2str(N) ".txt"]);
 %A1 = sum(A1);
 %A2 = sum(A2);
 %%AC = sum(AC)
%% end

%A = load("11.txt");
%B = load("22.txt");
%A = A(2:end);
%B = B(2:end);


%E = (A-B);
%max(abs(E))
%norm(E,2)/564
%max(abs(A-B))
%sum(A(2:end))
%4*pi
%sum(E)
%E
%sum(A)
%sum(B)
%F = load("F.txt");
%F = load("F1.txt");
%A = load("A1.txt");
%M = load("M1.txt");
%X = load("X1.txt");

%F = load("F2.txt");
%A = load("A2.txt");
%M = load("M2.txt");
%X = load("X2.txt");

%F = F(2:end);
%m = size(F,2);
%nume = size(A,1);
%X = X(2:end);
%area = sum(F);
%Aspr = sparse(A(:,1),A(:,2),A(:,3),m,m);
%Mspr = sparse(M(:,1),M(:,2),M(:,3),m,m);
%
%L = (Aspr+Mspr);
%res = norm((Aspr+Mspr)*X' - F');
%rAm = rank((Aspr+Mspr));
%cAm  = cond(Aspr+Mspr);
%x = (Aspr+Mspr)\F';
%res2 = norm((Aspr+Mspr)*x - F');
%gap = norm(X'-x);
%
%L - L'

F = load("rhsp2_iface.txt");
F = F(2:end);
m = size(F,2);
M = load("ap2_iface.txt");
Mspr = sparse(M(:,1),M(:,2),M(:,3),m,m);
area = sum((sum(Mspr)))


 
 
function [lambda,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,show,nrank,tau,sep,sepinf,DA,DB)

%SINGGEP   Finite eigenvalues of a singular matrix pencil
%
% [lambda,Z,d,X,Y,U,V,DA,DB] = singgep(A,B,show,tau,sep,sepinf,nrank,DA,DB)
% computes finite regular eigenvalues of a pencil (A,B) with a 
% rank-completing perturbation
%
% If matrices A, B are rectangular, we add zero columns or rows to make
% them square n x n matrices. Using one rank k random perturbation, where 
% k = n - nrank(A-lambda*B) we then compute all regular eigenvalues of 
% the pencil A-lambda*B
% 
% Input
%  - A,B : matrices of the same size n1 x n2
%  - show : display a table with values used for the separation or no (0)
%  - nrank : normal rank of A-lambda*B (if not given it is computed)
%  - tau : perturbation size (1e-2)
%  - sep : threshold for the separation of regular eigenvalues (sqrt(eps))
%  - sepinf : threshold for the separation of finite eigenvalues (100*eps)
%  - DA, DB : matrices of size k x k for k = max(n1,n2) - nrank for
%    prescribed eigenvalues (if not given random matrices are used)
% 
% Output
%  - lambda : finite regular eigenvalues found
%  - Z : matrix with criteria used for the separation (use show=1 to see it)
%         - col. 1: index 1 to n
%         - col. 2: eigenvalue of the perturbed pencil
%         - col. 3: |s(eigenvalue)|
%         - col. 4: ||V'*x||, where x is right eigenvector
%         - col. 5: ||U'*y||, where y is left eigenvector
%         - col. 6: max(||V'*x||, ||U'*y||)
%         - col. 7: regular eigenvalue (1) or not (0)
%         - col. 8: finite eigenvalue (1) or not (0)
%  - d : eigenvalues of the perturbed pencil
%  - X, Y: right and left eigenvectors of the perturbed pencil
%  - U, V: matrices with k orthonormal columns used for the perturbation
%  - DA, DB : matrices that define prescribed eigenvalues eig(DA,DB)

% Reference: Algorithm 1 in 
% M.E. Hochstenbach, C. Mehl, B. Plestenjak: Solving singular generalized 
% eigenvalue problems by a rank-completing perturbation, arXiv:1805:07657

% MultiParEig toolbox
% B. Plestenjak, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Bor Plestenjak
% 23.05.2018

% Matlab at least 2014a is required (we need left and right eigenvectors)
if verLessThan('matlab', '8.3')
    error('Matlab at least 2014a is required')
end

% multiple precision is supported
class_t = superiorfloat(A,B);

% Make sure all inputs are of the same numeric type.
if ~isa(A,class_t), A = numeric_t(A,class_t); end;
if ~isa(B,class_t), B = numeric_t(B,class_t); end;

[m,n] = size(A);

% if matrices are rectangular, we add zero columns or rows
if m~=n
    n = max(m,n);
    A(n,n) = 0;
    B(n,n) = 0;
end

% scaling of matrices so that ||A||_1=||B||_1=1
alfa = norm(A,1);
beta = norm(B,1);
A = A/alfa;
B = B/beta;

if nargin<3, show = 0; end
% if normal rank is not given, we compute it from a random linear combination of A and B 
if nargin<4, nrank = rank(rand(class_t)*A+rand(class_t)*B); end

k = n - nrank; 
U = orth(rand(n,k,class_t)); % random matrix with k orthonormal columns
V = orth(rand(n,k,class_t)); % random matrix with k orthonormal columns

if nargin<5, tau = numeric_t(1e-2,class_t); end
if nargin<6, sep = sqrt(eps(class_t)); end   % criteria for regular eigenvalues
if nargin<7, sepinf = 100*eps(class_t); end  % criteria for finite regular eigenvalues
% if matrices DA and DB are not given, we use diagonal matrices with
% elements randomly chosen from the interval [1,2]
if nargin<8
    DA = diag(numeric_t('1',class_t)+rand(k,1,class_t)); 
else
    if ~isa(DA,class_t), DA = numeric_t(DA,class_t); end;
end
if nargin<9
    DB = diag(numeric_t('1',class_t)+rand(k,1,class_t)); 
else
    if ~isa(DB,class_t), DB = numeric_t(DB,class_t); end;
end

% we apply rank k completing perturbation
A1 = A + tau*U*DA*V';
B1 = B + tau*U*DB*V';

% compute eigenvalues and left and right eigenvector (supported since Matlab 2014a)
[X,D,Y] = eig(A1,B1);

% eig does not return normalized eigenvectors, thus we normalize them
for i = 1:n
    X(:,i) = X(:,i)/norm(X(:,i));
    Y(:,i) = Y(:,i)/norm(Y(:,i));
end
d = diag(D)*alfa/beta; % eigenvalues (scaled back to original (A,B))
for i = 1:n
    s(i,1) = abs(Y(:,i)'*B1*X(:,i));   % 1/condition number
    a(i,1) = norm(V'*X(:,i));          % product V'*x
    b(i,1) = norm(U'*Y(:,i));          % product U'*y
    e(i,1) = max(a(i,1),b(i,1));       % criteria for separation
end

[esort,eord] = sort(e);  % we sort the eigenvalues by criteria for separation
d = d(eord);
s = s(eord);
a = a(eord);
b = b(eord);
X = X(:,eord);
Y = Y(:,eord);
e = esort;

pos1 = find(e<sep);     % subset of regular eigenvalues
[ssort,~] = sort(-s(pos1)); 
ssort = -ssort;

pos2 = find(ssort>sepinf); 

ind = intersect(find(e<sep),find(s>sepinf)); % subset of finite regular eigenvalues 
lambda = d(ind);

Z = [(1:n)' d s a b esort esort<sep s>sepinf]; 

if show
    disp(' ')
    disp(' i  |     re(lambda)     im(lambda) |     s(eig) |   ||V''*x|| |   ||U''*y|| |  criterion | reg. fin.')
    disp('---------------------------------------------------------------------------------------------------')
    for i=1:n
        j = i;
        fprintf('%3d | %14.6e %14.6e | %10.2e | %10.2e | %10.2e | %10.2e |  %d    %d\n',i,real(d(j)),imag(d(j)),s(j),a(j),b(j),e(j),e(j)<sep,s(j)>sepinf);
    end
end






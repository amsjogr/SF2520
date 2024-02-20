clc, clear all, close all
%%Part 1
%larger n, slower convergence
%Conj converges quicker, d=2, n=32,64,128
%d3, nvec 7,14 for theory for jacobi, quadruple. O(n^2). Conj double O(n)

% Part 1a
% K = 2000;
% dvec = [1,2,2]
% nvec = [32,64]
% 
% for u = dvec
% d = u;
% for i = nvec
% 
% 
% N = i^d;
% 
% % FIXED N
% % N = 4096;
% % i = round(10^(log10(N)/d))
% 
% A = lap(i,d);
% 
% b = rand(N,1);
% 
% 
% text2 = ['Conjugate, n = ', num2str(i),'d = ', num2str(u)];
% [x, res] = conj2(i, d, A, K, b);
% figure(1)
% semilogy(1:K, res,'DisplayName',text2,'LineWidth',2)
% hold on 
% 
% text1 = ['Jacobi, n = ', num2str(i), 'd = ', num2str(u)];
% [x2, res2] = jacob(A, i, d, K, b);
% semilogy(1:K-1, res2,'--','DisplayName',text1,'LineWidth',2)
% hold on
% 
% legend show 
% 
% end
% end


% Part 1b
% K = 2000;
% dvec = [2];
% nvec = [15,32,64];
% for u = dvec
% d = u;
% for i = nvec
% 
% A = lap(i,d);
% N = i^d;
% b = rand(N,1);
% disp('Conj')
% tic
% [x, res] = conj3(i, d, A, K, b);
% toc
% 
% disp('Backslash')
% tic
% xback =A\b;
% toc
% 
% end
% end


function [x, residuals] = conj3(n, d, A, K, b)
eps = 10 ^ -10;
N = n^d;
residuals=zeros(K,1);
x = zeros(N,1);
r = b;
p = x;
count = 1;

for i = 0:K-1
  if i == 0
       beta = 0;
  else
       beta = (r'*r)/(rold'*rold);
  end
       
  p = r + beta * p;
  pvec(:,count) = p;
  Ap = A*p;
  alpha = (p'*r)/(p'*Ap);
  x = x + alpha *p;
  rold = r;
  r = rold - alpha * Ap;
    
  res= norm((A*x-b)) / norm(b);
  residuals(count)=res;
  count = count + 1;
  if norm(r)/norm(b) < eps
      break
  end
    end
end 

function [x, residuals] = conj2(n, d, A, K, b)

N = n^d;
residuals=zeros(K,1);
x = zeros(N,1);
r = b;
p = x;
count = 1;

for i = 0:K-1
  if i == 0
       beta = 0;
  else
       beta = (r'*r)/(rold'*rold);
  end
       
  p = r + beta * p;
  pvec(:,count) = p;
  Ap = A*p;
  alpha = (p'*r)/(p'*Ap);
  x = x + alpha *p;
  rold = r;
  r = rold - alpha * Ap;
    
  res= norm((A*x-b)) / norm(b);
  residuals(count)=res;
  count = count + 1;
    end
end 

function [x res] = jacob(A, n, d, it, b)

M = diag(diag(A));
M = sparse(M);
T = M-A;
T = sparse(T);
N = n^d;
x = zeros(N,1);


for i=1:it-1

 x(:,i+1)=M\(T*x(:,i)+b);

 res(i)=norm((A*x(:,i+1)-b))/norm(b);

end


end


function A = lap(n,d)
% LAP
%    A = LAP(N,D) returns the system matrix corresponding to the 
%    standard second order discretization of the Laplace operator 
%    for D dimensions and N unknowns in each coordinate direction.
%    The size of the matrix is thus N^D times N^D.

e = ones(n,1);
A1 = -spdiags([e -2*e e], -1:1, n, n);
I1 = speye(n,n);

A=A1;
I=I1;
for k=2:d
  A = kron(A,I1)+kron(I,A1);
  I = kron(I,I1);
end
A = A*n^2;
end
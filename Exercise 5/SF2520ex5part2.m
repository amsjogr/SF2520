%2
clear all, close all, clc
%spy(A)


%% pcg part a & b
% load('cooling_flange.mat');
% tol = 10^-4;
% b = rand(length(A),1);
% maxit=2000;
% M=diag(diag(A));
% 
% tic
% [X,FLAG,RELRES,ITER,RESVEC] = pcg(A,b,tol,maxit);
% toc
% ITER
% RELRES
% figure(1)
% semilogy(1:length(RESVEC),RESVEC,'LineWidth',2)
% hold on 
% 
% tic
% [X,FLAG,RELRES,ITER,RESVEC] = pcg(A,b,tol,maxit, M);
% toc
% ITER
% RELRES
% semilogy(1:length(RESVEC),RESVEC,'LineWidth',2)
% hold on
% 
% LL = ichol(A);
% tic
% [X,FLAG,RELRES,ITER,RESVEC] = pcg(A,b,tol,maxit, LL, LL');
% toc
% ITER
% RELRES
% semilogy(1:length(RESVEC),RESVEC,'LineWidth',2)
% legend('pcg','pcg with M','pcg with ichol')
% tic
% xback =A\b;
% toc

%% gmres part c
% load('convdiff.mat');
% tol = 10^-4;
% b = rand(length(A),1);
% maxit=2000;
% 
% [X,FLAG,RELRES,ITER,RESVEC] = pcg(A,b,tol,maxit);
% FLAG %doesn't converge, FLAG not 0
% 
% [L,U] = ilu(A);
% M = diag(diag(A));
% 
% tic
% [X,FLAG,RELRES,ITER,RESVEC] = gmres(A,b,[],tol , maxit);
% toc
% ITER
% RELRES
% 
% figure(1)
% semilogy(1:length(RESVEC),RESVEC,'LineWidth',2)
% hold on
% 
% tic
% [X,FLAG,RELRES,ITER,RESVEC] = gmres(A,b,[],tol , maxit,M);
% toc
% ITER
% RELRES
% semilogy(1:length(RESVEC),RESVEC,'LineWidth',2)
% hold on
% 
% tic
% [X,FLAG,RELRES,ITER,RESVEC] = gmres(A,b,[],tol , maxit,L, U);
% toc
% ITER
% RELRES
% 
% semilogy(1:length(RESVEC),RESVEC,'LineWidth',2)
% 
% legend('gmres','gmres with M','gmres with ilu')
% 
% tic
% xback =A\b;
% toc
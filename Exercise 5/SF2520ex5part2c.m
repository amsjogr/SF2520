clc, clear all, close all
load convdiff.mat 
maxit = 2000;
tol = 10^-4;
b = rand(length(A),1);
[X,FLAG,RELRES,ITER,RESVEC] = pcg(A,b,tol,maxit);
FLAG %doesn't converge, FLAG not 0

% A_conv=A;
% a=A_conv'*A_conv;
% b=A_conv*A_conv'; %normality check
% 
% A_conv=sparse(A_conv);
% %spy(A_conv); 
% N=length(A_conv);
% b=rand(N,1);
% maxit=1500;
% tol=1e-4;
% 
% [X,FLAG,RELRES,ITER,RESVEC]=pcg(A_conv,b, tol, maxit);
% FLAGpcg= FLAG 
% maxit=500; %<<N=55096
% tic
% [X,FLAG,RELRES,ITER,RESVEC]=gmres(A_conv,b,[], tol, maxit); 
% gmrestime=toc;
% %experiments
% gmresFLAG=FLAG;gmresRELRES=RELRES,gmresITER=ITER(2),gmresRESVEC=RESVEC; %outer iteration
% figure()
% semilogy(RESVEC)
% title('RESVEC for gmres method')
% xlabel('iteration')
% ylabel('Residual size')
% 
% %Using backslash
% tic
% A_conv\b;
% 
% backslashtime_gmres=toc;
% 
% %b) PRECONDITION
% diagelements=diag(A_conv);
% M=diag(diagelements);
% tic
% [X,FLAG,RELRES,ITER,RESVEC]=gmres(A_conv,b,[], tol, maxit, M); %
% 
% precondtime_gmres=toc;
% %experiments
% precondFLAG=FLAG;precondRELRES=RELRES,precondITER=ITER(2),precondRESVEC=RESVEC;
% figure()
% semilogy(precondRESVEC)
% title('RESVEC for gmres method with preconditioner M=diag. part of A')
% xlabel('iteration')
% ylabel('Residual size')
% 
% 
% %Incomplete Cholesky
% [L, U]=ilu(A_conv);
% tic
% [X,FLAG,RELRES,ITER,RESVEC]=gmres(A_conv,b, [], tol, maxit, L, U); %
% ilutime_gmres=toc;
% %experiments
% iluFLAG=FLAG;iluRELRES=RELRES,iluITER=ITER(2),iluRESVEC=RESVEC;
% figure()
% semilogy(iluRESVEC)
% title('RESVEC for gmres method with incomplete LU factorization as preconditioner')
% xlabel('iteration')
% ylabel('Residual size')
% 
% gmrestime, backslashtime_gmres, precondtime_gmres, ilutime_gmres %print times
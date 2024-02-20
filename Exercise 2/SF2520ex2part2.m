clc, clear all, close all
Lx = 12;
Ly = 5;
Text = 25;
h = 0.2;

N = Lx/h+1;
M = Ly/h+1;

xvec = linspace(0,Lx,N);
yvec = linspace(0,Ly,M);

hx=zeros(N,1) + h;
hy=zeros(M,1) + h;


main=zeros(N,1); 
first_diag=zeros(N,1);
last_diag=zeros(N,1);

i=1:N;                                        
                           
 for j=1:M
     jj=i+(j-1)*N;
     main(jj)=(2./(2*hx')).*(-1./hx'+-1./hx')+(2/(2*h))*(-1/h+-1/h);    
     first_diag(jj)=(2./(2*hx')).*(1./hx');                                                          
     last_diag(jj)=(2/(2*h))*(1/h);                                                    
 end


T=zeros(N,M);  %Empty T matrix             
Tot=N*M;     

A=sparse(Tot,Tot);  %Empty A matrix
index=(1:Tot);
A=A+sparse(index,index,main,Tot,Tot); %Middle diagonal

diag_1=1:1:Tot-1;
diag_2=1:N:Tot-1;
diag_3=N:N:Tot-1;
diag_4=N-1:N:Tot-1;

A=A+sparse(diag_1,diag_1+1,first_diag(diag_1),Tot,Tot); 
A=A+sparse(diag_2,diag_2+1,first_diag(diag_2),Tot,Tot);
A=A+sparse(diag_3,diag_3+1,-first_diag(diag_3),Tot,Tot);
A=A+sparse(diag_1+1,diag_1,first_diag(diag_1+1),Tot,Tot); 
A=A+sparse(diag_4+1,diag_4,first_diag(diag_4+1),Tot,Tot);
A=A+sparse(diag_3+1,diag_3,-first_diag(diag_3+1),Tot,Tot);



diag5=N+1:1:Tot;
diag6=N+1:1:2*N;
A=A+sparse(diag5-N,diag5,last_diag(diag5-N),Tot,Tot); 
A=A+sparse(diag6-N,diag6,last_diag(diag6-N),Tot,Tot); 
A=A+sparse(diag5,diag5-N,last_diag(diag5),Tot,Tot); 
A=A+sparse(Tot-(diag6-N-1),Tot-diag6+1,last_diag(Tot-(diag6-N-1)),Tot,Tot); 

      
%f = zeros(N*M,1) - 2; %for f = 2

f=[]; 
    for i=1:M
        y=yvec(i);
        for k=1:N
            x=xvec(k);
            f = [f; heatsource(x,y)];
        end
    end

Nvec = [1:N]';
f(1:N)=Text; %Dirichlet boundary


A(Nvec,:)=0; 
A=A+sparse(Nvec,Nvec,1,Tot,Tot);

x=A\f; %Matrix solve

T=reshape(x,N,M);                
[xvec yvec]=meshgrid(xvec,yvec); %For plots           

%% Mesh plot
%mesh(xvec',yvec',T)
% xlabel('x')    
% ylabel('y')
% title('Temperatures of a 2D block')
%% Imagesc plot
%imagesc(T)
% camroll(90)
% xlabel('y')    
% ylabel('x')
% colorbar
%% Contour plot
% contour(T)
% camroll(90)
% xlabel('y')    
% ylabel('x')
% title('Contour plot')
% colorbar



function f = heatsource(x,y) %f function
f = -100 * exp(-0.5*(x-4).^2  - 4*(y-1).^2);
end
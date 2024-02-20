clc, clear all, close all
% vt + A*vx = F => vt = [ut;vt], vx = [ux,vx], F = [0,f]
%Hyperbolic: A = diagonalizable (independent eigenvalue vectors), real eigenvalues.
Fr = 0.35; %subcritical
%Fr = 1.5; %supercritical
alpha = 1/Fr;
A = [1, alpha;alpha, 1];
Hyp = eig(A);

N = 1000; %Grid points
L_0 = -0.4;
L_1 = 0.7;

dx = (L_1 - L_0)/N;
x = L_0 + [0:N]*dx;
%Stability check: lambda =< 1/max(abs(eig(A)))
c = 1/max(abs(eig(A))); %What is c?
dt = dx/max(abs(eig(A)));
% dt = dt * 1.01 %Stability check
%At N = 200, still stable at 1.01dt. Not stable at N = 1000.
%Lower N makes solution more lenient towards CFL
Tend = 0.15;
time = linspace(0,Tend,Tend/dt); %Time vector
lambda = dt/dx;

[heightslw, vellw] = LW(time,x,dx,dt,lambda,A,N);
[heightsup, velup] = UP(time,x,dx,dt,lambda,A,N);

figure(1) %Plotting u
hold on
plot(x,heightslw(end,:),'LineWidth',2)
plot(x,heightsup(end,:),'LineWidth',2,'LineStyle','--')
xlabel('x')
ylabel('u')
txt = ['Heights at t = ' num2str(Tend) ', N = ' num2str(N)];
title(txt)
legend('Lax-Wendroff','Upwind')

figure(2) %Plotting v
hold on
plot(x,vellw(end,:),'LineWidth',2)
plot(x,velup(end,:),'LineWidth',2,'LineStyle','--')
xlabel('x')
ylabel('velocity')
txt = ['Velocity at t = ' num2str(Tend) ', N = ' num2str(N)];
title(txt)
legend('Lax-Wendroff','Upwind')


%Peaks to the right are always same height, regardless of N or method
%Peaks to the left converge to same height with increased N. LW converge quicker. 

%Left wave length = 0.1 (ish)
%Right wave length = 0.2

%Upwind smears. LW oscillates. Especially breaking CFL.

%%3D plots
% figure(3) %3D plot u, Lax-Wendroff
% [xcords, tcords] = meshgrid(x,time);
% surf(xcords,tcords,heightslw(:,:))
% title('u approximated with Lax-Wendroff')
% xlabel('x'), ylabel('Time'), zlabel('u')
% figure(4)%3D plot u, Upwind
% surf(xcords,tcords,heightsup)
% xlabel('x'), ylabel('Time'), zlabel('u')
% title('u approximated with Upwind')

function [heightsLW, velocLW] = LW(time, x, dx, dt, lambda, A, N)
LW = zeros(2, N+1); %initial, row 1 = u height, row 2 = v velocity
heightsLW = LW(1,:);
velocLW = LW(2,:);
count = 2;
LWnew = LW; %initialize sol
f = @(x,t) (sin(40*pi*t+pi/6)>0.5).*(abs(x)<1/20).*sin(20*pi*x);
for t = time(1:end-1)
for i = 2:N - 1
F = [0; f( x(i),t ) + f( x(i),t+dt )] / 2 - (lambda/4) * A * [0; f( x(i+1),t ) - f( x(i-1)-dx,t )];
LWnew(:,i) = LW(:,i) - ((lambda/2) * A * ( LW(:,i+1) - LW(:,i-1) )) + (((lambda^2)/2) * A^2 * (LW(:,i+1) - 2 * LW(:,i) +LW(:,i-1))) + (dt * F); %inner points
end
LWnew(:,1) = 2 * LWnew(:,2) - LWnew(:,3);
LWnew(:,N) = 2 * LWnew(:, N-1) - LWnew(:,N-2); %boundaries interpolation
heightsLW(count,:) = LWnew(1,:);
velocLW(count,:) = LWnew(2,:);
count = count + 1;
LW = LWnew;
end
end
function [heightsUP, velocUP] = UP(time, x, dx, dt, lambda, A, N)
f = @(x,t) (sin(40*pi*t+pi/6)>0.5).*(abs(x)<1/20).*sin(20*pi*x);
[S L]=eig(A); 
Lp = L.*(L>0); 
Lm = L.*(L<0);
Am=S*Lm*inv(S); 
Ap=S*Lp*inv(S);
upnew = zeros(2,N+1); %initialize sol
up = upnew;
count = 2;
for t = time(1:end-1)
for i = 2:N - 1
F = [0; f(x(i),t)];
upnew(:,i) = up(:,i) - (lambda * Ap * ( up(:,i) - up(:,i-1) )) - lambda* Am * (up(:,i+1) - up(:,i)) + dt * F; %inner points
end
upnew(:,1) = 2 * upnew(:,2) - upnew(:,3);
upnew(:,N) = 2 * upnew(:, N-1) - upnew(:,N-2); %boundaries interpolation
heightsUP(count,:) = upnew(1,:);
velocUP(count,:) = upnew(2,:);
count = count + 1;
up = upnew;
end
end

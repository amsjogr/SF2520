clear all, close all, clc

v = 1;
N = [10, 20, 40, 80];
h = ones(1,length(N))./N;


%% a)

for i = h
[z,u] = FDM(i,v);
plot(z,u)
hold on 
%Temps(i) = u(N(i)/2 + 1);
end 
xlabel('z')
ylabel('Temperature')
legend('N = 10', 'N = 20', 'N = 40', 'N = 80')

%% b)
%v = [1, 5, 15, 100];
%for ii = 1:length(N)
 %   h = 1/N(ii);
% for i = 1:length(v)
% [z,T] = FDM(h,v(i));
% subplot(2,2,ii)
% pp = ['N = ' num2str(N(ii))];
% title(pp)
% set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',16);
% xlabel('z')
% ylabel('Temperature')
% hold on 
% end 
% legend('v = 1', 'v = 5', 'v = 15', 'v = 100')
%end


%%error
% error = [];
% v = 1;
% for i = h
% [z, u] = solve(i,v);
% [zcorr,ucorr] = solve(i/2,v);
% error = [error, abs(u(end)-ucorr(end))];
% end 
% 
% loglog(h, error)
% hold on
% loglog(h, h.^2)
% legend('Error', 'h^2')
% title('Error on log log scale')
% xlabel('h')





function [z, u] = FDM(h,v)

T_out = 25;
T_0 = 100;
L = 1;
a = 0.1;
b = 0.4;
Q_0 = 7000;
alpha_0 = 50;



alpha = @(v) sqrt(((v^2)/4) + (alpha_0^2)) - (v/2);
z = linspace(h,1,1/h);
ind = length(z);


 adiag = -1/(h^2) - (v/(2*h));
 bdiag = 2/(h^2);
 cdiag = -1/(h^2) + (v/(2*h));


A = diag([bdiag+zeros(1,ind-1) 0], 0);
A = A + diag(cdiag+zeros(1,ind-1), 1);
A = A + diag([adiag+zeros(1,ind-2) 0], -1);
A(ind, ind-1) = adiag+cdiag;
A(ind, ind) = bdiag - cdiag*(2*h*alpha(v));      



f = [];
for i = z
if i >= a && i <= b
        f = [f Q_0*sin((i-a)*pi/(b-a))];
else 
    f = [f 0];
end
end
f(1) = f(1) - adiag*T_0;
f(end) = f(end) - cdiag * 2 * h * alpha(v) * T_out;
u = A\f';
u = [T_0; u];
z = [0 z]';
    
end


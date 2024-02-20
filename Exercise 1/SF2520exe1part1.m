%SF2520 exercise 1
%Part 1
clc, clear all, close all

N = [50, 100, 200, 400, 800, 1600];

alpha = 0.07;
a = (1/4) * [1,sqrt(11),2]';
T = 50;
m_0 = [0,0,1]';

cc = 1;

for n = N
t_vector = linspace(0,T,n);
h = T/n;
m = [];
m(:,1) = m_0;
count = 1;
hh(cc) = h;
for i = t_vector

   m_n = m(:,count);
   m(:,count + 1) = RK(m_n,h,a,alpha);

    count = count + 1;

end

norm = vecnorm(m).^(-1);
mfin = m*diag(norm);
hold on
plot3([mfin(1,:)], [mfin(2,:)], [mfin(3,:)])
results(:,cc) = m(:,end) ;

cc = cc + 1;

end

%diff = vecnorm(results(:,1:5)-results(1,2:6));
%loglog(diff)

% figure(1)
% subplot(3,1,1);
% plot(t_vector,mfin(1,2:end))
% title('X-coordinates over time')
% xlabel('Time')
% ylabel('X position')
% 
% subplot(3,1,2); 
% plot(t_vector,mfin(2,2:end))
% title('Y-coordinates over time')
% xlabel('Time')
% ylabel('Y position')
% 
% subplot(3,1,3); 
% plot(t_vector,mfin(3,2:end))
% title('Z-coordinates over time')
% xlabel('Time')
% ylabel('Z position')
% 
% 
% figure(2)
% plot3([0,a(1)],[0,a(2)],[0,a(3)])
% axis equal
% hold on 
%plot3([mfin(1,:)], [mfin(2,:)], [mfin(3,:)])
% legend('a-vector','norm m pos')


function [m_new] = RK(m_n,h,a,alpha)

%k1
k_1 = cross(a, m_n) + alpha .* (cross(a, cross(a, m_n)));

%k2
m_2 = m_n + h.*k_1;
k_2 = cross(a, m_2) + alpha .* (cross(a, cross(a, m_2)));

%k3
m_3 = m_n + (h/4)*k_1 + (h/4).*k_2;
k_3 = cross(a, m_3) + alpha .* (cross(a, cross(a, m_3)));

m_new = m_n + (h/6) .* (k_1 + k_2 + 4 .* k_3);

end



%SF2520 exercise 1
%Part 2
clc, clear all, close all

my = 1/82.45;
B = [0 1;-1 0];
r_0 = [-my, 0]'; %earth
r_1 = [1 - my, 0]'; %moon
r = [1.15,0]';
rprim =  [0, -0.975]';


r_5 = [0.4681;0.6355];
r_20 = [-0.2186;-0.2136];
r_40 = [-1.4926; -0.3339];


% Tend = 20;
% h = Tend/24500;
% 
% [r,rprim] = RK(my,r,r_0,r_1,h,rprim,B);
% [r,rprim] = AdamsB(r, rprim, r_0, r_1, B, my, Tend, h);
% figure(1)
% plot(r(1,:), r(2,:))
% hold on
% axis equal
% scatter(r_0(1),r_0(2),"filled")
% scatter(r_1(1),r_1(2),"filled")
% 
% 
% legend('h = 0.005', 'h = 0.001','h = 0.0005','earth','moon')
% title('Trajectory 2D w/ different h')
% xlabel('x coords')
% ylabel('y coords')

% figure(1)
% plot3(r(1,:), r(2,:), t_vec)
% title('Position coordinates')
% xlabel('x')
% ylabel('y ')
% zlabel('Time')
% legend('Trajectory')
% 
% figure(2)
% plot3(rprim(1,:), rprim(2,:), t_vec)
% title('Velocities')
% xlabel('x')
% ylabel('y')
% zlabel('Time')
% legend('r prim')


options = odeset('RelTol', 0.001);
initials = [1.15,0,0,-0.975]';
[t, y] = ode23(@sat4ode, [0 20], initials, options);


% plot(y(:,1),(y(:,2)))
% scatter(r_0(1),r_0(2),"filled")
% scatter(r_1(1),r_1(2),"filled")
% legend('AdamsB4', 'Ode23','earth','moon')
% title('AdamsB4 and Ode23 trajectory')
% xlabel('x coords')
% ylabel('y coords')

% stepdiff = diff(t);
% maxstep = max(stepdiff)
% minstep = min(stepdiff)
% figure(2)
% plot(t(1:end-1), stepdiff)
% title('Step size over time')
% ylabel('Step size')
% xlabel('Time')




function dydt = sat4ode(t,r)
my = 1/82.45;
B = [0 1;-1 0];
r_0 = [-my, 0]'; 
r_1 = [1 - my, 0]'; 
gg = (-(1-my)*(([r(1); r(2)]-r_0)./(vecnorm([r(1); r(2)]-r_0).^3))) - (my*(([r(1); r(2)]-r_1)./(vecnorm([r(1); r(2)]-r_1).^3))) + (2*B*[r(3); r(4)]) + [r(1); r(2)];
dydt = [r(3), r(4), gg(1), gg(2)]';
end 

function [r,rprim] = RK(my,r,r_0,r_1,h,rprim,B)

sat = @(r,rprim) [rprim,(-(1-my)*((r-r_0)./(vecnorm(r-r_0).^3))) - (my*((r-r_1)./(vecnorm(r(:,1)-r_1).^3))) + (2*B*rprim) + r];

for i = 1:3

k_1 = sat(r(:,i), rprim(:,i));
k_2 = sat(r(:,i) + h * k_1(:,1), rprim(:,i) + h * k_1(:,2));
k_3 = sat(r(:,i) + (h/4)*k_1(:,1) +  (h/4)*k_2(:,1), rprim(:,i) + (h/4)*k_1(:,2) + (h/4)*k_2(:,2));


rprim(:,i+1) = rprim(:,i) + (h/6) * (k_1(:,2)+k_2(:,2) + 4 * k_3(:,2));
r(:,i+1) = r(:,i) + (h/6) * (k_1(:,1) + k_2(:,1) + 4 * k_3(:,1));

end
end


function [r, rprim] = AdamsB(r, rprim, r_0, r_1, B, my, Tend, h)

for i = 4:Tend/h-1
sat = @(r,rprim) [rprim,(-(1-my)*((r-r_0)./(vecnorm(r-r_0).^3))) - (my*((r-r_1)./(vecnorm(r(:,1)-r_1).^3))) + (2*B*rprim) + r];


    f_1 = sat(r(:,i), rprim(:,i));

    f_2 = sat(r(:,i-1), rprim(:,i-1));
    f_3 = sat(r(:,i-2), rprim(:,i-2));
    f_4 = sat(r(:,i-3), rprim(:,i-3));

    r(:,i+1) = r(:,i) + (h/24)*(55*f_1(:,1) - 59*f_2(:,1) + 37*f_3(:,1) - 9*f_4(:,1));
    rprim(:,i+1) = rprim(:,i) + (h/24)*(55*f_1(:,2) - 59*f_2(:,2) + 37*f_3(:,2) - 9*f_4(:,2));
end
end

function [r, rprim] = eulerex(r, rprim, tend, h, r_0, r_1, B, my)
sat = @(r,rprim) [rprim,(-(1-my)*((r-r_0)./(vecnorm(r-r_0).^3))) - (my*((r-r_1)./(vecnorm(r(:,1)-r_1).^3))) + (2*B*rprim) + r];

for i = 1:tend/h
k1 = sat(r(:,i), rprim(:,i));
r(:,i+1) = r(:,i) + h*k1(:,1);
rprim(:,i+1) = rprim(:,i) + h*k1(:,2);
end
end


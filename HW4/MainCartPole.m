%{
Zane Grothe
AERO 7350
HW 4
11/10/23
%}

clear all
close all
clc

%{
This script (accompanied with function file TPBVP_2ndOrder.m) finds the
minimum energy dynamics for a cart-pole swing-up and balance for a given
amount of time.
Part 1 derives the generalized necessary conditions for optimality.
Part 2 solves the 2pt BVP for the specified boundary conditions.
This file creates 3 function files and returns 3 figures for analytical
comparison.
%}


%% Part 1

% Constants
m1 = 1;   % kg
m2 = 0.3; % kg
g = 9.81; % m/s^2
l = 0.5;  % meters
tf = 1; % time duration for minimum energy trajectory

NofEq = 4; % establish size of our problem
syms t x th xd thd u

% Create states and costates vectors
xbar = [x;th;xd;thd]; % states
lambar = cell(NofEq,1);
for i = 1:NofEq
    lambar{i} = sprintf('lam%d',i);
end
lambar = lambar(:);
lambar = sym(lambar,'real'); % costates

L = 1/2*u^2; % define the Lagrangian

% Equations of motion
xdd = (m2*l*thd^2*sin(th) + m2*g*cos(th)*sin(th) + u) ...
                            / (m1 + m2*(1-cos(th)^2));
thdd = -(m2*l*thd^2*cos(th)*sin(th) + (m1+m2)*g*sin(th) + u*cos(th)) ...
                                       / (m1*l + m2*l*(1-cos(th)^2));
% State dynamics
xbard = [xd;thd;xdd;thdd];

%% Hamiltonian
H = L + lambar.'*xbard;

%% Costate Dynamics
lambard = -jacobian(H,xbar).';

%% Optimal Control
u0 = solve(diff(H,u),u);

w = [xbard;lambard];
w = subs(w,u,u0);
z = [xbar,lambar];

H0 = subs(H,u,u0);

% Create files for calculating these values to be called later
fc = matlabFunction(w,'file','CartPoleEOM.m','vars',{t,z},'outputs',{'zdot'});
dc = matlabFunction(u0,'file','Find_u0.m','vars',{z},'outputs',{'u0'});
ec = matlabFunction(H0,'file','Find_H0.m','vars',{z},'outputs',{'H0'});

%% Part 2

% Boundary Conditions
x0 = zeros(4,1); % start hanging at rest
xf = zeros(4,1); % finish balancing at rest
xf(2) = pi;

tspan = [0,tf];

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
options_fsolve = optimoptions('fsolve','Display','iter',...
    'MaxFunEvals',5000,'MaxIter',1000,'TolFun',1e-10,'Tolx',1e-10);
% use fsolve to iterate for a solution to the costate values
norm_fval = 1;
while norm_fval > 1e-4
    lambar0 = rand(4,1); % initial random guess
    [lam_sol,fval,exitflag,output] = fsolve(@TPBVP_2ndOrder,lambar0,options_fsolve,x0,xf,tspan,options);
    norm_fval = norm(fval); % compare the norms of the vectors
end

%% Plots

% Establish the final solution vector and propogate optimal trajectory
z0 = [x0;lam_sol];
[t_final,z_final] = ode45(@CartPoleEOM,tspan,z0,options);
Nd = length(t_final);
u0 = zeros(1,Nd);
H0 = u0;

% Calculate control and Hamiltonian for established final solution
for i = 1:Nd
    z_current = z_final(i,:).';
    u0(i) = Find_u0(z_current);
    H0(i) = Find_H0(z_current);
end

% Extract states and costates
x = z_final(:,1);
th = z_final(:,2);
xd = z_final(:,3);
thd = z_final(:,4);
lam1 = z_final(:,5);
lam2 = z_final(:,6);
lam3 = z_final(:,7);
lam4 = z_final(:,8);

figure(1) % States

subplot(2,2,1);
plot(t_final,x)
title('Cart Position')
ylabel('x (m)')

subplot(2,2,2);
plot(t_final,th)
title('Pendulum Angle')
ylabel('theta (rad)')

subplot(2,2,3);
plot(t_final,xd)
title('Cart Velocity')
ylabel('xdot (m/s)')

subplot(2,2,4);
plot(t_final,thd)
title('Pendulum Angular Velocity')
ylabel('theta dot (rad/s)')

figure(2) % Control and Hamiltonian

subplot(2,1,1);
plot(t_final,u0)
title('Control')
ylabel('u (Newtons)')

subplot(2,1,2);
plot(t_final,H0)
title('Hamiltonian')

figure(3) % Costates

subplot(2,2,1);
plot(t_final,lam1)
title('Costate 1')

subplot(2,2,2);
plot(t_final,lam2)
title('Costate 2')

subplot(2,2,3);
plot(t_final,lam3)
title('Costate 3')

subplot(2,2,4);
plot(t_final,lam4)
title('Costate 4')


%% Animation

% Constants
a = 30; % cart width (cm)
b = 15; % cart height (cm)
d = 5; % cart wheels diameter (cm)
ground = 1; % ground thickness (cm)

% Convert previous linear dimensions from meters to centimeters (*100)
x = x * 100; % cart position (cm)
l = l * 100; % pendulum length (cm)

% Adjust pendulum angle by pi/2
theta = th - pi/2;

% Cart position
Ax = x;
Ay = 0;
AM = max(abs(Ax));

time = length(t_final);
axislimitsfactor = 1.5;
alf = axislimitsfactor;
axislimits = [-AM*alf AM*alf -l*alf l*alf]; % keep the cart and pendulum on plot
%%
for w = 1:time
    figure(4)
    groundX = [-(AM+a/2)*1.1,(AM+a/2)*1.1,(AM+a/2)*1.1,-(AM+a/2)*1.1];
    groundY = [Ay-b/2-d-ground,Ay-b/2-d-ground,Ay-b/2-d,Ay-b/2-d];
    fill(groundX,groundY,'b') % ground polygon (blue)
    hold on
    cartX = [Ax(w)-a/2 Ax(w)+a/2 Ax(w)+a/2 Ax(w)-a/2];
    cartY = [Ay-b/2 Ay-b/2 Ay+b/2 Ay+b/2];
    fill(cartX,cartY, 'r') % cart polygon (red)
    pos1 = [Ax(w)-a/2+d Ay-b/2-d d d]; % position for wheel 1
    rectangle('Position',pos1,'Curvature',[1,1]) % wheel 1 polygon
    pos2 = [Ax(w)+a/2-2*d Ay-b/2-d d d]; % position for wheel 2
    rectangle('Position',pos2,'Curvature',[1,1]) % wheel 2 polygon
    plot([Ax(w),l*cos(theta(w))+Ax(w)],[Ay,l*sin(theta(w))+Ay],'k','linewidth',2) % pendulum
    p = plot(l*cos(theta(w))+Ax(w),l*sin(theta(w))+Ay,'ko'); % pendulum bob
    p.MarkerFaceColor = 'g';
    p.MarkerSize = 8/l*100;
    p.MarkerEdgeColor = 'g';
    axis(axislimits) % so the axis doesn't change
    axis equal
    axis off % keep the axis dimensions but make them hidden
    hold off
end

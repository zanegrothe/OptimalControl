%{
Zane Grothe
AERO 7350
Final Project
12/5/23


This script (accompanied with function file MI_TransferTraj.m) uses the 
variational approach (indirect method) to transfer a rocket vehicle from
a given initial circular orbit to the largest possible circular orbit for
a given amount of time.
Part 1 derives the generalized necessary conditions for optimality.
Part 2 solves the 2pt BVP for the specified boundary conditions.
Part 3 plots the solutions
This file creates 3 function files and returns 3 figures for analytical
comparison.
%}

clear
close all
clc

%% Part 1 Necessary Conditions

% Constants
t0       = 0; % Initial time (days)
tf       = 193; % Time of flight (days)
m0       = 4535.9; % Initial spacecraft mass (kg)
T        = 3.781; % Thrust (N)
mdot     = 5.85; % Propellant flow rate (kg/day)
mu_Sun   = 132712441933; % Sun's gravitational constant (km^3/s^2)
r0_Earth = 149597870.691; % Radius of Earth's circular orbit (km)
a_Sun    = mu_Sun/r0_Earth^2; % Gravitational acceleration of the Sun on Earth (km/s^2)
a_p      = T/m0/1e3; % Acceleration produced by the propulsion system at the departure from Earth (km/s^2)

aN       = a_p/a_Sun; % Non-dimentional acceleration due to the thrusting
tfN      = tf/sqrt(r0_Earth^3/mu_Sun)*24*60*60; % Non-dimensional total flight time
muN      = 1; % Non-dimentional gravitational constant
mdotN    = mdot/m0*tf/tfN; % Non-dimensional propellant flow rate

% Establish size of our problem
NofEq = 3;
syms t r u v phi th

% Create states and costates vectors
xbar = [r;u;v]; % States
lambar = cell(NofEq,1);
for i = 1:NofEq
    lambar{i} = sprintf('lam%d',i);
end
lambar = lambar(:);
lambar = sym(lambar,'real'); % Costates

L = 0; % Define the Lagrangian

Ta = aN/(1-mdotN*t); % Augmented thrust

% Equations of motion
rdot = u;
udot = v^2/r - muN/r^2 + Ta*sin(phi);
vdot = -u*v/r + Ta*cos(phi);
thetadot = v/r;

% State dynamics
xbardot = [rdot;udot;vdot];

% Hamiltonian
H = L + lambar.'*xbardot;

% Costate Dynamics
lambardot = -jacobian(H,xbar).';

% Optimal Control
%phi0 = min(solve(diff(H,phi),phi)); Does not solve quadrant ambiguity
phi0 = atan2(-lambar(2),-lambar(3));
sinphi = -lambar(2)/sqrt(lambar(2)^2+lambar(3)^2);
cosphi = -lambar(3)/sqrt(lambar(2)^2+lambar(3)^2);

% Substitute optimal control into state and costate dynamics
y = [xbardot;lambardot];
y = subs(y,sin(phi),sinphi);
y = subs(y,cos(phi),cosphi);
z = [xbar;lambar];

% Substitute optimal control into state and costate dynamics (including theta)
ytheta = [xbardot;lambardot;thetadot];
ytheta = subs(ytheta,sin(phi),sinphi);
ytheta = subs(ytheta,cos(phi),cosphi);
ztheta = [xbar;lambar;th];

% Substitute optimal control into the Hamiltonian
H0 = subs(H,sin(phi),sinphi);
H0 = subs(H0,cos(phi),cosphi);

% Create files for calculating these values to be called later
dc = matlabFunction(y,'file','RocketEOM.m','vars',{t,z},'outputs',{'zdot'});
fc = matlabFunction(H0,'file','Find_H0.m','vars',{t,z},'outputs',{'H0'});
gc = matlabFunction(ytheta,'file','RocketEOMandT','vars',{t,ztheta},'outputs',{'zdottheta'});

%% Part 2 TPBVP

% Boundary Conditions
r0N = 1;  % Non-dimensional initial radial distance
u0N = 0;  % Non-dimensional initial radial velocity
v0N = 1;  % Non-dimensional initial transverse velocity
th0N = 0; % Non-dimensional initial radial angle
x0N  = [r0N;u0N;v0N]; % Initial BC vector

ufN = 0;    % Non-dimensional final radial velocity
xfN  = ufN; % Final BC vecor using transversality

tspan = [t0,tfN]; % Time duration for maximum radial distance trajectory

options = odeset('RelTol',1e-8,'AbsTol',1e-8);
options_fsolve = optimoptions('fsolve','Display','iter',...
    'MaxFunEvals',1000,'MaxIter',1000,'TolFun',1e-10,'Tolx',1e-10);

% Use fsolve to iterate for a solution to the costate values
norm_fval = 1;
tic
while norm_fval > 1e-4
    lambar0 = rand(3,1); % Initial random guess
    [lam_sol,fval,exitflag,output] = fsolve(@MI_TransferTraj,lambar0,options_fsolve,x0N,xfN,tspan,options);
    norm_fval = norm(fval); % Compare the norms of the vectors
end
toc

%% Part 3 Plots

% Establish the final solution vector and propogate optimal trajectory
z0 = [x0N;lam_sol;th0N];
[t_final,z_final] = ode45(@RocketEOMandT,tspan,z0,options);

% Extraction
r = z_final(:,1);
u = z_final(:,2);
v = z_final(:,3);
phi = atan2(-z_final(:,5),-z_final(:,6));
neg_phi = phi < 0;
phi(neg_phi) = phi(neg_phi) + 2*pi;
theta = z_final(:,7);

% Figures
% Plot r,u,v over time
figure(1)
subplot(3,1,1)
plot(t_final,r)
xlabel('Time')
ylabel('r')
title('Radial Distance')
subplot(3,1,2)
plot(t_final,u)
xlabel('Time')
ylabel('u')
title('Radial Velocity')
subplot(3,1,3)
plot(t_final,v)
xlabel('Time')
ylabel('v')
title('Transverse Velocity')

% Plot control over time
figure(2)
plot(t_final(1:end-1),phi(1:end-1))
xlabel('Time')
ylabel('phi')
title('Thrust Steering Angle')

% Plot final trajectory solution
figure(3)
z = linspace(-r(1),r(1));
yi = sqrt(r(1)^2-z.^2);
plot(z,yi,'b--')
hold on
w = linspace(-r(end),r(end));
yf = sqrt(r(end)^2-w.^2);
plot(w,yf,'g--')
trajx = r.*cos(theta);
trajy = r.*sin(theta);
plot(trajx,trajy,'k-')
plot(trajx(1),trajy(1),'or')
plot(trajx(end),trajy(end),'xr')
xlim([-r(end)-0.5,r(end)+0.5])
ylim([0,max(yf)+0.1])
xlabel('X (AU)')
ylabel('Y (AU)')
title('Final Orbit Transfer Trajectory')
legend('Initial Orbit','Final Orbit','Orbit Transfer','Departure','Injection')


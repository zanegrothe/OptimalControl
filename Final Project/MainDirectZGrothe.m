%{
Zane Grothe
AERO 7350
Final Project
12/5/23


This script (accompanied with function file MD_Dynf.m) uses the
optimization toolbox in MATLAB to solve a double-integrator problem using
the problem-based approach. The transcription method used is
Hermite-Simpson. The problem solved is the transfer of a rocket vehicle 
from a given initial circular orbit to the largest possible circular orbit 
for a given amount of time.
Part 1 transcribes the problem to a discrete set of points
Part 2 solves the double-integrator problem using the optimization function fmincon
Part 3 plots the solutions
This file returns 3 figures for analytical comparison.
%}

clear
close all
clc

%% Part 1 Transcription

% Transcription method
Transmethod = 'HermiteSimpson';

% Boundary conditions
Xi = [1;0;1;0]; % [r,u,v,theta]

% Define decision variables (states, controls and any other unknown parameter)
N = 50; % number of mesh/grid points
X = optimvar('X',4,N,'LowerBound',0,'UpperBound',10);
U = optimvar('U',1,N,'LowerBound',-pi,'UpperBound',pi);
rf = optimvar('rf',1,'LowerBound',1,'UpperBound',30);

% If we wanted, we could define X and give a different name to its
% components and also define separate lower and upper limits for each
% X = [optimvar('x',1,N,'LowerBound',-10,'UpperBound',10);
%      optimvar('xdot',1,N,'LowerBound',-20,'UpperBound',20)];

% Define the time interval
tvec = linspace(0,3.32,N);

% Define the optmization problem: give it a name and determine that we want
% to minimize the objective
doubleintprob = optimproblem('ObjectiveSense','minimize');

% Construct defect equations
deffectscons = optimconstr(4,N-1);

switch Transmethod
    case 'HermiteSimpson'
        for i=2:N
            tk   = tvec(i-1);
            tkp1 = tvec(i);
            hs   = tkp1-tk;
            Xk   = X(:,i-1);
            Xkp1 = X(:,i);
            Uk   = U(:,i-1);
            Ukp1 = U(:,i);
            fk   = MD_Dynf(tk,Xk,Uk);
            fkp1 = MD_Dynf(tkp1,Xkp1,Ukp1);
            tm   = 0.5*(tk+tkp1);
            Xm   = 0.5*(Xk+Xkp1)+hs/8*(fk-fkp1);
            Um   = 0.5*(Uk+Ukp1);
            fm   = MD_Dynf(tm,Xm,Um);
            deffectscons(:,i) = Xkp1-Xk-1/6*hs*(fk+4*fm+fkp1)==0; % note the "==" to define equality constraint
        end
end

% Apply a state-only inequality constraint (based on the solution of the unconstrained solution)
% stateineqcon = optimconstr(1,N);
% for i=1:N
   % stateineqcon(i) = X(1,i) >= 1;  % here we set a lower bound on the first state
% end
% doubleintprob.Constraints.stateineq = stateineqcon; % add the constraint to the field of constraints

% Apply a control-only inequality constraint (based on the solution of the unconstrained solution)
controlineqcon = optimconstr(1,N);
for i=1:N
   controlineqcon(i) = U(i) <= pi;  % Here we set an upper bound on the control
end
doubleintprob.Constraints.controlineq = controlineqcon; % Add the constraint to the field of constraints

% Add defects to the problem "Constraints" field
doubleintprob.Constraints.deffects = deffectscons;

% Form boundary constraints
State_i_con = X(:,1)   == Xi; 
State_f_con1 = X(2,end) == 0;
State_f_con2 = X(3,end) == sqrt(1/X(1,end));

% Add boundary conditions to the problem "Constraints" field
doubleintprob.Constraints.State_i_con = State_i_con;
doubleintprob.Constraints.State_f_con1 = State_f_con1;
doubleintprob.Constraints.State_f_con2 = State_f_con2;

% Define objecive
costproblem = -rf;
doubleintprob.Objective = costproblem; % add to "Objective" field

% Initial guess
Xg = repmat(Xi,1,N);
Ug = repmat(-1,1,N);
rfg = 2;
init.X = Xg;
init.U = Ug;
init.rf = rfg;

%% Part 2 Double-Integrator Problem

% solve (it uses Automatic Differentiation, read the problem-based details)
% solve only provides first-order derivatives, but we can augment the
% solver with Hessian information. Again, you can read on help.

options = optimoptions('fmincon','Algorithm','interior-point',...
                       'TolX',1e-5,...
                       'TolFun',1e-4,...
                       'TolCon',1e-5,...
                       'Display','iter');
tic                   
[sol_ds,fval,exitflag,output] = solve(doubleintprob,init,'options', options,ConstraintDerivative="finite-differences",ObjectiveDerivative="finite-differences");
toc
beep

%% Part 3 Plots

% Extraction
States = sol_ds.X;
Controls = sol_ds.U;

r     = States(1,:);
u     = States(2,:);
v     = States(3,:);
theta = States(4,:);

neg_U = Controls < 0;
Controls(neg_U) = Controls(neg_U) + 2*pi;

% Figures
% Plot r,u,v over time
figure(4)
subplot(3,1,1)
plot(tvec,r)
xlabel('Time')
ylabel('r')
title('Radial Distance')
subplot(3,1,2)
plot(tvec,u)
xlabel('Time')
ylabel('u')
title('Radial Velocity')
subplot(3,1,3)
plot(tvec,v)
xlabel('Time')
ylabel('v')
title('Transverse Velocity')

% Plot control over time
figure(5)
plot(tvec(1:end-1),Controls(1:end-1))
xlabel('Time')
ylabel('phi')
title('Thrust Steering Angle')

% Plot final trajectory solution
figure(6)
z = linspace(-r(1,1),r(1,1));
yi = sqrt(r(1,1)^2-z.^2);
plot(z,yi,'b--')
hold on
w = linspace(-r(1,end),r(1,end));
yf = sqrt(r(1,end)^2-w.^2);
plot(w,yf,'g--')
trajx = r.*cos(theta);
trajy = r.*sin(theta);
plot(trajx,trajy,'k-')
plot(trajx(1),trajy(1),'or')
plot(trajx(end),trajy(end),'xr')
axis fill
xlim([-r(1,end)-0.5,r(1,end)+0.5])
ylim([0,max(yf)+0.1])
xlabel('X (AU)')
ylabel('Y (AU)')
title('Final Orbit Transfer Trajectory')
legend('Initial Orbit','Final Orbit','Orbit Transfer','Departure','Injection')


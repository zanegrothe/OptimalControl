%{
Author: Ehsan Taheri
Department of Aerospace Engineering
Auburn University
Last modified: 11/24/2020

This code uses MATLAB's problem-based approach to solve the minimum-time double-integrator problem.
Problem-based approach is much more intuitive as you define the constraint
and objective as we write them mathematically. The "solver" will
automatically determine the arguments needed for solving the problem

Transcription Method: Euler-Forward 
 
Problem: A minimum-time problem for a second-order system. 

minimize J = t_f

subject to 

State dynamics:     xddot = u
Control bounds:     -1<= u <=1
Initial conditions:  x(0) = 10, xdot(0) = 0 
Final conditions:    x(tf) = 0, xdot(tf) = 0

%}

clear; clc

% transcription method
Transmethod = 'EulerForward';

% define boundary conditions
Xi = [10;0];
Xf = [0;0];

% define decision variables (states, controls and any other unknown parameter)
N = 60; % number of mesh/grid points
X = optimvar('X',2,N);
U = optimvar('U',1,N,'LowerBound',-1,'UpperBound',1);
tf = optimvar('tf',1,'LowerBound',1,'UpperBound',10);

% if we wanted, we could define X and give a different name to its
% components and also define separate lower and upper limits for each
% X = [optimvar('x',1,N,'LowerBound',-10,'UpperBound',10);
%      optimvar('xdot',1,N,'LowerBound',-20,'UpperBound',20)];

% define the time interval
tau = linspace(0,1,N);
tvec = tau*tf;

% define the optmization problem: give it a name and determine that we want
% to minimize the objective
doubleintprob = optimproblem('ObjectiveSense', 'minimize');

% construct defect equations
deffectscons = optimconstr(2,N-1);

switch Transmethod
    case 'EulerForward'
        for i=2:N
            hs = tvec(i)-tvec(i-1);
            Xkp1 = X(:,i);
            Xk   = X(:,i-1);
            Uk   = U(:,i-1);
            tk   = tvec(i-1);
            deffectscons(:,i) = Xkp1-Xk-hs*Dynf(tk,Xk,Uk)==0; % note the "==" to define equality constraint
        end
end

% add defects to the problem "Constraints" field
doubleintprob.Constraints.deffects = deffectscons;

% form boundary constraints
State_i_con = X(:,1)   == Xi; 
State_f_con = X(:,end) == Xf;

% add boundary conditions to the problem "Constraints" field
doubleintprob.Constraints.State_i_con = State_i_con;
doubleintprob.Constraints.State_f_con = State_f_con;

%% define objecive
costproblem = tf; 
doubleintprob.Objective = costproblem; % add to "Objective" field

%% initial guess
Xg = repmat(Xi,1,N);
Ug = repmat(0,1,N);
tfg = 6;
init.X = Xg;
init.U = Ug;
init.tf = tfg;

%% solve (it uses Automatic Differentiation, read the problem-based details)
% solve only provides first-order derivatives, but we can augment the
% solver with Hessian information. Again, you can read on help.

options = optimoptions('fmincon','Algorithm','interior-point',...
                       'TolX',1e-5,...
                       'TolFun',1e-4,...
                       'TolCon',1e-5,...
                       'Display','iter');

tic                   
[sol_ds,fval,exitflag,output] = solve(doubleintprob,init,'options', options);
toc

%% Plot section
t_vec = tau*sol_ds.tf;
States = sol_ds.X;
Controls = sol_ds.U;

x1 = States(1,:);
x2 = States(2,:);

figure(3)
hold on
plot(t_vec,x1,'.-')
plot(t_vec,x2,'.-')
grid off
legend('$x_1$','$x_2$','interpreter','latex')
xlabel('Time')

figure(4)
hold on
plot(t_vec(1:end-1),Controls(1:end-1),'o-')
grid off
xlabel('Time')
ylabel('Control')

figure(5)
plot(States(1,:),States(2,:),'.-')
hold on
plot(States(1,1),States(2,1),'go')
plot(States(1,end),States(2,end),'ro')
grid 
xmax = max(States(1,:));

xlim([-xmax,xmax])
ylim([-xmax,xmax])
axis square 
xlabel('x1')
ylabel('x2')
legend('Solution','Initial','Final')
%{
Zane Grothe
AERO 7350
Final Project
12/5/23

   MI_TransferTraj.m is written to form a residual vector to be minimized 
by fsolve in MainIndirectZGrothe.m. This file is used as part of a single 
shooting method that propagates the states and costates forward using the
equations of motion in RocketEOM.m through the time span given. The
ending values (zf) are then compared to the constraints of the problem
(Gamf) and the file passes the residual vector.

Input:
   lambar0  : Initial guess for costate vector [3x1] (lambda bar initial)
Constants:
   x0N      : Initial conditions for the state vector [3x1] (x initial normalized)
   xfN      : Final condition for the state vector [1x1] (x final normalized)
   tspan    : Initial and final time for ode45 [1x2] (time span)
   options  : Relative and absolute tolerances for ode45
Ouput:
   Residual : Vector 
%}

function Residual = MI_TransferTraj(lambar0,x0N,xfN,tspan,options)

z0 = [x0N;lambar0];

[~,z_guess] = ode45(@RocketEOM,tspan,z0,options);
zf = z_guess(end,:).';

x_guess = [zf(4);...
           zf(2);...
           zf(3)];
Gamf = [zf(6)/2*sqrt(1/zf(1))^3-1;...
                              xfN;...
                    sqrt(1/zf(1))];

Residual = x_guess - Gamf;

end

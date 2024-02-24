function Xdot = MD_Dynf(t,X,U) 
%{
State dynamics file
Inputs
    t = current time
    X = current state vector
    U = current control vector
Outputs
   xdot = rhs equations of motion
%}

Ta = 0.1405/(1-0.07497*t); % Augmented Thrust

rdot     = X(2);
udot     = X(3)^2/X(1) - 1/X(1)^2 + Ta*sin(U);
vdot     = -X(2)*X(3)/X(1) + Ta*cos(U);
thetadot = X(3)/X(1);

Xdot = [rdot;udot;vdot;thetadot];
end
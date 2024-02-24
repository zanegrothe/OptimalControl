function Xdot = Dynf(~,X,U) 
%{
State dynamics file
Inputs
    t = current time
    X = current state vector
    U = current control vector
Outputs
   xdot = rhs equations of motion
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x        = X(1); % position
xdot     = X(2); % velocity 

u = U(1); % control
Xdot = [xdot;
        u];
end
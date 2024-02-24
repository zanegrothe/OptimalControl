function Residual = TPBVP_2ndOrder(lambar0,x0,xf,tspan,options)

z0 = [x0;lambar0];

[~,z_guess] = ode45(@CartPoleEOM,tspan,z0,options);
zf = z_guess(end,:).';

x_guess = zf(1:4);

Residual = x_guess - xf;

end

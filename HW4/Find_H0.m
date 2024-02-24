function H0 = Find_H0(in1)
%Find_H0
%    H0 = Find_H0(IN1)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    08-Feb-2024 16:39:00

lam1 = in1(5);
lam2 = in1(6);
lam3 = in1(7);
lam4 = in1(8);
th = in1(2);
thd = in1(4);
xd = in1(3);
t2 = cos(th);
t3 = sin(th);
t4 = thd.^2;
t6 = lam3.*1.0e+1;
t5 = t2.^2;
t8 = lam4.*t2.*2.0e+1;
t7 = t5.*3.0;
t9 = -t8;
t10 = t7-1.3e+1;
t11 = t6+t9;
t12 = 1.0./t10;
H0 = (t11.^2.*t12.^2)./2.0+lam2.*thd+lam1.*xd+(lam4.*(t3.*1.2753e+1+t2.*t3.*t4.*(3.0./2.0e+1)+t2.*t11.*t12))./(t5.*(3.0./2.0e+1)-1.3e+1./2.0e+1)-(lam3.*(t2.*t3.*2.943+t3.*t4.*(3.0./2.0e+1)+t11.*t12))./(t5.*(3.0./1.0e+1)-1.3e+1./1.0e+1);
end

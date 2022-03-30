function xdot = duffingDE(X,omega,beta)
%% States and Costates
x1 = X(1);
x2 = X(2);
L1 = X(3);
L2 = X(4);
%% Dynamic equations
x1dot = x2;
x2dot = -(omega^2*x1 + beta*x1.^3 + L2);
L1dot = L2.*(omega^2+3*beta*x1.^2);
L2dot = - L1;
xdot = [x1dot; x2dot; L1dot; L2dot];
end



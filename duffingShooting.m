function [R, J] = duffingShooting(lambdas,omega,beta,BC,t,BCtype)
x10 = BC(1);
x20 = BC(2);
x1f = BC(3);
x2f = BC(4);
L1 = lambdas(1);
L2 = lambdas(2);
ICs_shooting = [x10,x20,L1,L2]; 
fshooting = @(t,x) duffingDE(x,omega,beta);
opts = odeset('RelTol',1e-13,'AbsTol',1e-13);
[T,X] = ode45(fshooting,t,ICs_shooting,opts);

% switch BCtype
%     case 'fixed'
%         R = [X(end,1) - x1f; X(end,2) - x2f];
%     case 'free'
%         R = [L1(end) - (X(end,1) - x1f); L2(end) - (X(end,2) - x2f)];
%         J = eye(2);
%     case 'P0-Pf'
%         R = [X(end,1) - x1f; L2(end) - (X(end,2)-x2f)];
% end

        R = [L1(end) - (X(end,1) - x1f); L2(end) - (X(end,2) - x2f)];

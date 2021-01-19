function J = duffingJacobian_ipopt(X,BC,omega,beta,D,N)
%% States and Costates
x1 = X(1:N);
x2 = X(N+1:2*N);
L1 = X(2*N+1:3*N);
L2 = X(3*N+1:4*N);
x10 = BC(1);
x20 = BC(2);
x1f = BC(3);
x2f = BC(4);

%% Residual
R(1:N) = D*x1 - x2;
R(N+1:2*N) = D*x2 + omega^2*x1 + beta*x1.^3 + L2;
R(2*N+1:3*N) = D*L1 - L2.*(omega^2 + 3*beta*x1.^2);
R(3*N+1:4*N) = D*L2 + L1;
% BCs
R(N+1)   = x1(1) - x10;
R(2*N+1) = x2(1) - x20;
R(2*N)   = x1(end) - x1f;
R(3*N)   = x2(end) - x2f;




%% Jacobian
dR1dx1 = D; %
dR1dx2 = -ones(N);
dR1dL1 = zeros(N);
dR1dL2 = zeros(N);
%         dR2dx1 = omega^2+2*beta*x1*x1';
dR2dx1 = omega^2 + 3*beta*diag(x1.^2);
dR2dx2 = D;
dR2dL1 = zeros(N);
dR2dL2 = ones(N);
dR3dx1 = -6*beta*diag(L2.*x1);
dR3dx2 = zeros(N);
dR3dL1 = D;
%         dR3dL2 = -(omega^2+3*beta*x1*x1');
dR3dL2 = -(omega^2 + 3*beta*diag(x1.^2));

dR4dx1 = zeros(N);
dR4dx2 = zeros(N);
dR4dL1 = ones(N);
dR4dL2 = D;



J = [dR1dx1 dR1dx2 dR1dL1 dR1dL2;...
    dR2dx1 dR2dx2 dR2dL1 dR2dL2;...
    dR3dx1 dR3dx2 dR3dL1 dR3dL2;...
    dR4dx1 dR4dx2 dR4dL1 dR4dL2];


J(N+1,:) = [ones(1,N), zeros(1,3*N)];
J(2*N+1,:) = [zeros(1,N), ones(1,N), zeros(1,2*N)];
J(2*N,:) = [ones(1,N), zeros(1,3*N)];
J(3*N,:) = [zeros(1,N), ones(1,N), zeros(1,2*N)];

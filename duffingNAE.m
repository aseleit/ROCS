function [R,J] = duffingNAE(X,BC,omega,beta,D,N,BCtype)
%% States and Costates
x1 = X(1:N);
x2 = X(N+1:2*N);
L1 = X(2*N+1:3*N);
L2 = X(3*N+1:4*N);
x10 = BC(1);
x20 = BC(2);
x1f = BC(3);
x2f = BC(4);
s11 = 1;
s22 = 1;

% Residual
R(1:N) = D*x1 - x2 ;
R(N+1:2*N) = D*x2 + omega^2*x1 + beta*x1.^3 + L2;
R(2*N+1:3*N) = D*L1 - L2.*(omega^2 + 3*beta*x1.^2);
R(3*N+1:4*N) = D*L2 + L1;

switch BCtype
    
    case "fixed"
        % BCs
        R(N+1)   = x1(1) - x10;
        R(2*N+1) = x2(1) - x20;
        R(2*N)   = x1(end) - x1f;
        R(3*N)   = x2(end) - x2f;
        % Jacobian
        dR1dx1 = D; dR1dx2 = -eye(N); dR1dL1 = zeros(N); dR1dL2 = zeros(N);
        dR2dx1 = diag(omega^2+3*beta*x1.^2); dR2dx2 = D; dR2dL1 = zeros(N); dR2dL2 = eye(N);
        dR3dx1 = diag(-6*beta*x1.*L2); dR3dx2 = zeros(N); dR3dL1 = D; dR3dL2 = diag(-omega^2-3*beta*x1.^2);
        dR4dx1 = zeros(N); dR4dx2 = zeros(N); dR4dL1 = eye(N); dR4dL2 = D;
        % Constraints
        dR2dx1(1,:) = [1,zeros(1,N-1)]; dR2dx2(1,:) = zeros(1,N);       dR2dL1(1,:) = zeros(1,N); dR2dL2(1,:) = zeros(1,N);
        dR2dx1(N,:) = [zeros(1,N-1),1]; dR2dx2(N,:) = zeros(1,N);       dR2dL1(N,:) = zeros(1,N); dR2dL2(N,:) = zeros(1,N);
        dR3dx1(1,:) = zeros(1,N);       dR3dx2(1,:) = [1,zeros(1,N-1)]; dR3dL1(1,:) = zeros(1,N); dR3dL2(1,:) = zeros(1,N);
        dR3dx1(N,:) = zeros(1,N);       dR3dx2(N,:) = [zeros(1,N-1),1]; dR3dL1(N,:) = zeros(1,N); dR3dL2(N,:) = zeros(1,N);
        
    case "free"
        %BCs
        R(1)   = x1(1) - x10;
        R(N+1) = x2(1) - x20;
        R(3*N) = L1(end) - s11*(x1(end)-x1f);
        R(4*N) = L2(end) - s22*(x2(end)-x2f);
        % Jacobian
        dR1dx1 = D; dR1dx2 = -eye(N); dR1dL1 = zeros(N); dR1dL2 = zeros(N);
        dR2dx1 = diag(omega^2+3*beta*x1.^2); dR2dx2 = D; dR2dL1 = zeros(N); dR2dL2 = eye(N);
        dR3dx1 = diag(-6*beta*x1.*L2); dR3dx2 = zeros(N); dR3dL1 = D; dR3dL2 = diag(-omega^2-3*beta*x1.^2);
        dR4dx1 = zeros(N); dR4dx2 = zeros(N); dR4dL1 = eye(N); dR4dL2 = D;
        % Constraints
        dR1dx1(1,:) = [1,zeros(1,N-1)];  dR1dx2(1,:) = zeros(1,N);        dR1dL1(1,:) = zeros(1,N);       dR1dL2(1,:) = zeros(1,N);
        dR2dx1(1,:) = zeros(1,N);        dR2dx2(1,:) = [1,zeros(1,N-1)];  dR2dL1(1,:) = zeros(1,N);       dR2dL2(1,:) = zeros(1,N);
        dR3dx1(N,:) = [zeros(1,N-1),-1]; dR3dx2(N,:) = zeros(1,N);        dR3dL1(N,:) = [zeros(1,N-1),1]; dR3dL2(N,:) = zeros(1,N);
        dR4dx1(N,:) = zeros(1,N);        dR4dx2(N,:) = [zeros(1,N-1),-1]; dR4dL1(N,:) = zeros(1,N);       dR4dL2(N,:) = [zeros(1,N-1),1];
        
    case "P0-Pf"
        %BCs %working
%                 R(1)   = x1(1) - x10;
%                 R(N+1) = L2(1);
%                 R(3*N) = x1(end) - x1f;
%                 R(4*N) = L2(end) - s22*(x2(end)-x2f);
%         
        R(1)   = x1(1) - x10;
        R(2*N+1) = L2(1);
        R(N) = x1(end) - x1f;
        R(3*N+1) = L2(end) - s22*(x2(end)-x2f);
        % Jacobian
        dR1dx1 = D; dR1dx2 = -eye(N); dR1dL1 = zeros(N); dR1dL2 = zeros(N);
        dR2dx1 = diag(omega^2+3*beta*x1.^2); dR2dx2 = D; dR2dL1 = zeros(N); dR2dL2 = eye(N);
        dR3dx1 = diag(-6*beta*x1.*L2); dR3dx2 = zeros(N); dR3dL1 = D; dR3dL2 = diag(-omega^2-3*beta*x1.^2);
        dR4dx1 = zeros(N); dR4dx2 = zeros(N); dR4dL1 = eye(N); dR4dL2 = D;
        % The Constraints
%         dR1dx1(1,:) = [1,zeros(1,N-1)]; dR1dx2(1,:) = zeros(1,N);        dR1dL1(1,:) = zeros(1,N);        dR1dL2(1,:) = zeros(1,N);
%         dR2dx1(N,:) = zeros(1,N);       dR2dx2(N,:) = zeros(1,N);        dR2dL1(N,:) = zeros(1,N);        dR2dL2(N,:) = [1,zeros(1,N-1)];
%         dR3dx1(N,:) = [zeros(1,N-1),1];       dR3dx2(1,:) = zeros(1,N);        dR3dL1(1,:) = zeros(1,N);        dR3dL2(1,:) = zeros(1,N);
%         dR4dx1(N,:) = zeros(1,N);       dR3dx2(N,:) = [zeros(1,N-1),-1]; dR3dL1(N,:) = zeros(1,N);        dR3dL2(N,:) = [zeros(1,N-1),1]; 
        
%         % The Constraints
        dR1dx1(1,:) = [1,zeros(1,N-1)]; dR1dx2(1,:) = zeros(1,N);        dR1dL1(1,:) = zeros(1,N);        dR1dL2(1,:) = zeros(1,N);
        dR1dx1(N,:) = [zeros(1,N-1),1]; dR1dx2(N,:) = zeros(1,N);        dR1dL1(N,:) = zeros(1,N);        dR1dL2(N,:) = zeros(1,N);
        dR3dx1(1,:) = zeros(1,N);       dR3dx2(1,:) = zeros(1,N);        dR3dL1(1,:) = zeros(1,N);        dR3dL2(1,:) = [1, zeros(1,N-1)];
        dR3dx1(N,:) = zeros(1,N);       dR3dx2(N,:) = [zeros(1,N-1),-1]; dR3dL1(N,:) = zeros(1,N);        dR3dL2(N,:) = [zeros(1,N-1),1];
        
    case "P0-Vf"
        %BCs
        R(1)   = x1(1) - x10;
        R(N+1) = L2(1);
        R(3*N) = L1(end) - s11*(x1(end)-x1f);
        R(4*N) = x2(end) - x2f;
    case "V0-Pf"
        %BCs
        R(1)   = L1(1);
        R(N+1) = x2(1) - x20;
        R(3*N) = x1(end) - x1f;
        R(4*N) = L2(end) - s22*(x2(end)-x2f);
        %
        
end


J = [dR1dx1, dR1dx2, dR1dL1, dR1dL2; ...
    dR2dx1, dR2dx2, dR2dL1, dR2dL2; ...
    dR3dx1, dR3dx2, dR3dL1, dR3dL2;...
    dR4dx1, dR4dx2, dR4dL1, dR4dL2];



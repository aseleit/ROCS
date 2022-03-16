% Solution of the Nonlinear duffing oscillator problem 
% The optimal control problem is approximated using CRBFS
% The system of NAEs is solved using: 
% 1. fsolve (working solution) - minimization of residual vector
% 2. IPOPT (not working so far)
% 3. fmincon - minimization of a scalar
format long
clear all; close all; clc;
epsilon = 1;
N = 40;% Number of nodes
M = 1; % Number of elements
Neq = 4;
tf = 2;
omega = 1; %rad/sec
x10 = 0;
x20 = 0;
x1f = 2;
x2f = 1;
BCtype = 'fixed';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'
BC = [x10,x20,x1f,x2f]';
beta = 0.9;
% initial_guess = 1:1:16;
% initial_guess = [2.35*ones(1*N*M,1); 1.3*ones(1*N*M,1); 2.35*ones(1*N*M,1); 1.3*ones(1*N*M,1)]; 
initial_guess = 0*ones( 4*N*M,1);
% initial_guess = [1*ones(N*M,1);3.4*ones(N*M,1);1*ones(N*M,1);3.8*ones(N*M,1)];    %works for global collocation
% initial_guess = [1*rand(2*N*M,1);.8*rand(2*N*M,1)];
%         initial_guess = [2.35*ones(2*N*M,1);1.3*ones(2*N*M,1)];
%% Approximation matrix
D = zeros(N,N,M);
phi = zeros(N,N,M);
phid = zeros(N,N,M);
for k = 1 : M
    t0e            = (k-1)*tf/M; % Initial time of each segment (element)
    tfe            = k*tf/M;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    %     tinterp(:,:,k) = linspace(t0e,tfe,5*n);
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
        %         phi_interp(i,:,k) = rbf0(c,te(:,i,k),tinterp(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end

%% IPOPT Solution
% ip=1;
% duffingJacobian_ipopt(initial_guess,BC,omega,beta,D,N); %to build the Jac Sparsity str
% 
% funcs.objective = @(X) duffingNAE_local(X,BC,omega,beta,D,N,M,ip);
% funcs.gradient = @(X) duffingGradient_ipopt(X,BC,omega,beta,D,N); %the gradient of the objective function
% 
% funcs.iterfunc  = @callback;
% 
% ipopt_options.lb = -2*ones(1,4*N);
% ipopt_options.ub = [5*ones(1,N), 2*ones(1,N), 5*ones(1,2*N)];
% ipopt_options.cl = zeros(1,4*N);             % Lower bounds on constraints.
% ipopt_options.cu = zeros(1,4*N);             % Upper bounds on constraints.
% 
% funcs.constraints = @(XX) constraints(XX,BC,omega,beta,D,N,M,ip);
% funcs.jacobian = @(X) duffingJacobian_ipopt(X,BC,omega,beta,D,N); % The Jacobian (1st derivative) of the Constraints
% funcs.jacobianstructure = @()duffingJacobianStr_ipopt(N);
% 
% ipopt_options.ipopt.hessian_approximation      = 'limited-memory';%limited-memory
% ipopt_options.ipopt.limited_memory_update_type = 'bfgs';
% ipopt_options.ipopt.derivative_test            = 'first-order';
% 
% 
% 
% % Set the IPOPT ipopt_options.
% ipopt_options.ipopt.linear_solver =  'mumps'; %'ma57'; 'pardiso';
% ipopt_options.ipopt.mu_strategy   = 'adaptive';
% ipopt_options.ipopt.print_level   = 5;
% ipopt_options.ipopt.tol           = 1e-14;
% ipopt_options.ipopt.max_iter      = 1000;
% 
% % Run IPOPT.
% [x_ipopt, info] = ipopt(initial_guess,funcs,ipopt_options);
%% CRBF Solution using fmincon
% fmincon_options = optimset('Display','iter','TolX',1e-14);
% [x_fmincon,Fval,ExitFlag,Output] = fmincon(@(x) duffingNAE_local(x,BC,omega,beta,D,N,M,ip),initial_guess,[],[],...
%     [],[],[],[],@(XX) constraints(XX,BC,omega,beta,D,N,M,ip),fmincon_options); %if constraints exist..
% [x_fmincon,Fval,ExitFlag,Output] = fmincon(@(x) duffingNAE_local(x,BC,omega,beta,D,N,M,ip),initial_guess,[],[],...
%     [],[],[],[],[],fmincon_options); %if constraints don't exist..

%% CRBF Solution using fsolve
ip=0;

Jlocal_pattern = ones(4*N);
J_pattern = kron(eye(M),Jlocal_pattern);

J_local = zeros(Neq*N,Neq*N,M);
J    = zeros(Neq*N*M, Neq*N*M); % Local matrices are to be concatenated sequentially
ns = 0;
Elb   = zeros(1,M);
Erb   = Elb;
SLb   = zeros(Neq,1,M);
SRb   = SLb;
R = zeros(N,M);
sizeJ_local = size(J_local);
LJ_local    = sizeJ_local(1);


%% Newton's method
[resid,JacB] = duffingNAE_local(initial_guess,BC,omega,beta,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
solu = initial_guess;

naes = 'NEWTON'; naep = 0.3;
tol = 1e-9; MaxIter = 100; iter = 0; istop = 0;
compTime_CRBFnewton_start = tic;
while(istop==0)
    if(iter>MaxIter)
        disp('Max iteration reached for the nonliner solver')
        istop = 1;
    end
    dsolu = nae_su(resid,JacB,solu,initial_guess,naes,naep,iter); % compute incremental solution   
    solu = solu + dsolu; 
    [resid,JacB] = duffingNAE_local(solu,BC,omega,beta,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
    % determine convergence
    if(norm(resid)<tol)
        istop = 1;
    end
    iter = iter + 1;
    fprintf('Iteration = %d\n',iter)
end
compTime_CRBFnewton_end = toc(compTime_CRBFnewton_start);
SOLU = reshape(solu,[4*N,M]);
X1_newton = SOLU(1:N,:);
X2_newton = SOLU(N+1:2*N,:);
L1_newton = SOLU(2*N+1:3*N,:);
L2_newton = SOLU(3*N+1:4*N,:);

X1_newton = X1_newton(:);
X2_newton = X2_newton(:);
L1_newton = L1_newton(:);
L2_newton = L2_newton(:);

ucrbf_newton = -L2_newton;
J_newton = 0.5*trapz(ucrbf_newton.^2);




f = @(X) duffingNAE_local(X,BC,omega,beta,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
options = optimoptions(@fsolve,'Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',true,'StepTolerance',1e-20,'FunctionTolerance',1e-20,'UseParallel',true,'FiniteDifferenceType','forward');
% options = optimoptions(@fsolve,'Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',false,'StepTolerance',1e-20,'FunctionTolerance',1e-20,'UseParallel',true,'JacobPattern',J_pattern,'FiniteDifferenceType','central');

[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);
% [xx,fval1,exitflag1,output1] = lsqnonlin(f,initial_guess);
x = reshape(xx,[4*N,M]);
X1 = x(1:N,:);
X2 = x(N+1:2*N,:);
L1 = x(2*N+1:3*N,:);
L2 = x(3*N+1:4*N,:);


%% Removing repeated values from T and states
T = te(:);
[~, ind] = unique(T);
duplicate_ind = flip(setdiff(1:size(T), ind));
duplicate_value = T(duplicate_ind);
for i = 1 : length(duplicate_ind)
    T(duplicate_ind(i))  = [];
    X1(duplicate_ind(i)) = [];
    X2(duplicate_ind(i)) = [];
    L1(duplicate_ind(i)) = [];
    L2(duplicate_ind(i)) = [];
end
%% Validation
% IC = [X1(1) X2(2)    3.401297538141612    0.659726264262221];
IC = [X1(1) X2(1) L1(1) L2(1)];

fexact = @(t,Xexact)duffingDE(Xexact,omega,beta);
opts = odeset('RelTol',1e-20,'AbsTol',1e-20);
[t, Xexact] = ode23(fexact,T,IC,opts);
X1error = abs(X1' - Xexact(:,1));
% X1error = abs(X1 - x1iclocs);
X2error = abs(X2' - Xexact(:,2));
% X2error = abs(X2 - x2iclocs);
L1error = abs(L1' - Xexact(:,3));
L2error = abs(L2' - Xexact(:,4));
Errors = [X1error,X2error];
% Uerror = abs(L2 - uiclocs);
%% Plotting
% figure(1)
plot(t,X1,'*-')
hold on
plot(t,X2,'+-')
% plot(t,-L2,'*')
axis('tight')
% legend('x1','x2','u')
grid on
plot(t,Xexact(:,1),'o')
hold on
plot(t,Xexact(:,2),'o')
% plot(t,-Xexact(:,4),'o')
axis('tight')
legend('x1','x2','u','x1_exact','x2_exact','u_exact')
xlabel('time[sec]')
ylabel('Exact states and control')
for pp = 1:M-1
    tline = te(:,end,pp);
    xline(tline,'linewidth',2);
end
hold off
figure(2)
plot(T,L1)
hold on
plot(T,L2)
grid on
legend('lambda_1','lambda_2')
for pp = 1:M-1
    tline = te(:,end,pp);
    xline(tline,'linewidth',2);
end
hold off

%%
function b = callback(t, f, x)
  fprintf('%3d  %0.3g \n',t,f);
  b = true;
end


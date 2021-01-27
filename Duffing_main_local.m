% Solution of the Nonlinear duffing oscillator problem 
% The optimal control problem is approximated using CRBFS
% The system of NAEs is solved using: 
% 1. fsolve (working solution) - minimization of residual vector
% 2. IPOPT (not working so far)
% 3. fmincon - minimization of a scalar

clear all; close all; clc;
epsilon = 1;
N = 40;% Number of nodes
M = 1; % Number of elements
tf = 2;
omega = 1; %rad/sec
x10 = 0;
x20 = 0;
x1f = 5;
x2f = 2;
% BCtype = 'fixed';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'
BC = [x10,x20,x1f,x2f]';
beta = 0.9;
% initial_guess = [1*ones(1*N,1);3.5*ones(1*N,1);3.5*ones(1*N,1);3.5*ones(1*N,1)]; 
% initial_guess = 1*ones(4*N*M,1);
initial_guess = [1*ones(2*N*M,1);3.8*ones(2*N*M,1)];
% initial_guess = [1*rand(2*N*M,1);.8*rand(2*N*M,1)];

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
        D(:,:,k) = phid(:,:,k)/phi(:,:,k);

end

%% IPOPT Solution
ip=1;
funcs.objective = @(X) duffingNAE_local(X,BC,omega,beta,D,N,M,ip);
funcs.gradient = @(X) duffingGradient_ipopt(X,BC,omega,beta,D,N); %the gradient of the objective function

funcs.iterfunc  = @callback;

%  ipopt_options.cl = [ 0 0 0 0 ];             % Lower bounds on constraints.
%  ipopt_options.cu = [ 5 2 inf inf];             % Upper bounds on constraints.

% funcs.constraints = @(XX) constraints(XX,BC,omega,beta,D,N,M,ip);
% funcs.jacobian = @(X) duffingJacobian_ipopt(X,BC,omega,beta,D,N); % The Jacobian (1st derivative) of the Constraints
% funcs.jacobianstructure = @()duffingJacobianStr_ipopt(N);

ipopt_options.ipopt.hessian_approximation      = 'limited-memory';%limited-memory
ipopt_options.ipopt.limited_memory_update_type = 'bfgs';
ipopt_options.ipopt.derivative_test            = 'first-order';



% Set the IPOPT ipopt_options.
ipopt_options.ipopt.linear_solver =  'mumps'; %'ma57'; 'pardiso';
ipopt_options.ipopt.mu_strategy   = 'adaptive';
ipopt_options.ipopt.print_level   = 5;
ipopt_options.ipopt.tol           = 1e-14;
% ipopt_options.ipopt.max_iter      = 1000;

% Run IPOPT.
[x_ipopt, info] = ipopt(initial_guess,funcs,ipopt_options);

%% CRBF Solution using fmincon
fmincon_options = optimset('Display','iter','TolX',1e-14);
% [x_fmincon,Fval,ExitFlag,Output] = fmincon(@(x) duffingNAE_local(x,BC,omega,beta,D,N,M,ip),initial_guess,[],[],...
%     [],[],[],[],@(XX) constraints(XX,BC,omega,beta,D,N,M,ip),fmincon_options); %if constraints exist..
[x_fmincon,Fval,ExitFlag,Output] = fmincon(@(x) duffingNAE_local(x,BC,omega,beta,D,N,M,ip),initial_guess,[],[],...
    [],[],[],[],[],fmincon_options); %if constraints don't exist..

%% CRBF Solution using fsolve
ip=0;
f = @(X) duffingNAE_local2(X,BC,omega,beta,D,N,M,ip);
options = optimset('Display','iter','Algorithm','Levenberg-Marquardt','Jacobian','off','TolX',1e-20,'TolFun',1e-20);
[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);

X1 = xx(1:N*M);
X2 = xx(N*M+1 : 2*N*M);
L1 = xx(2*N*M+1 : 3*N*M);
L2 = xx(3*N*M+1 : 4*N*M);

% X1 = [xx(1:N); xx(4*N+1:5*N)];
% X2 = [xx(N+1:2*N); xx(5*N+1:6*N)];
% L1 = [xx(2*N+1:3*N); xx(6*N+1:7*N)];
% L2 = [xx(3*N+1:4*N); xx(7*N+1:8*N)];


%% Removing repeated values from T and states
T = te(:);
[~, ind] = unique(T);
duplicate_ind = flip(setdiff(1:size(T), ind));
duplicate_value = T(duplicate_ind);
for i = 1 : length(duplicate_ind)
    T(duplicate_ind(i)) = [];
    X1(duplicate_ind(i)) = [];
    X2(duplicate_ind(i)) = [];
    L1(duplicate_ind(i)) = [];
    L2(duplicate_ind(i)) = [];
end

figure(1)
plot(T,X1)

figure(2)
plot(T,X2)
%%
function b = callback(t, f, x)
  fprintf('%3d  %0.3g \n',t,f);
  b = true;
end


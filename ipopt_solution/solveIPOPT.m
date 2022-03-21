function solverIPOPT()


% IPOPT Solution
ip=1;
duffingJacobian_ipopt(initial_guess,BC,omega,beta,D,N); %to build the Jac Sparsity str

funcs.objective = @(X) duffingNAE_local(X,BC,omega,beta,D,N,M,ip);
funcs.gradient = @(X) duffingGradient_ipopt(X,BC,omega,beta,D,N); %the gradient of the objective function

funcs.iterfunc  = @callback;

ipopt_options.lb = -2*ones(1,4*N);
ipopt_options.ub = [5*ones(1,N), 2*ones(1,N), 5*ones(1,2*N)];
ipopt_options.cl = zeros(1,4*N);             % Lower bounds on constraints.
ipopt_options.cu = zeros(1,4*N);             % Upper bounds on constraints.

funcs.constraints = @(XX) constraints(XX,BC,omega,beta,D,N,M,ip);
funcs.jacobian = @(X) duffingJacobian_ipopt(X,BC,omega,beta,D,N); % The Jacobian (1st derivative) of the Constraints
funcs.jacobianstructure = @()duffingJacobianStr_ipopt(N);

ipopt_options.ipopt.hessian_approximation      = 'limited-memory';%limited-memory
ipopt_options.ipopt.limited_memory_update_type = 'bfgs';
ipopt_options.ipopt.derivative_test            = 'first-order';



% Set the IPOPT ipopt_options.
ipopt_options.ipopt.linear_solver =  'mumps'; %'ma57'; 'pardiso';
ipopt_options.ipopt.mu_strategy   = 'adaptive';
ipopt_options.ipopt.print_level   = 5;
ipopt_options.ipopt.tol           = 1e-14;
ipopt_options.ipopt.max_iter      = 1000;

% Run IPOPT.
[x_ipopt, info] = ipopt(initial_guess,funcs,ipopt_options);
% CRBF Solution using fmincon
fmincon_options = optimset('Display','iter','TolX',1e-14);
[x_fmincon,Fval,ExitFlag,Output] = fmincon(@(x) duffingNAE_local(x,BC,omega,beta,D,N,M,ip),initial_guess,[],[],...
    [],[],[],[],@(XX) constraints(XX,BC,omega,beta,D,N,M,ip),fmincon_options); %if constraints exist..
[x_fmincon,Fval,ExitFlag,Output] = fmincon(@(x) duffingNAE_local(x,BC,omega,beta,D,N,M,ip),initial_guess,[],[],...
    [],[],[],[],[],fmincon_options); %if constraints don't exist..

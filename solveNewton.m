% ***************************************************************************************************
% Function solverNewton() {this function could be a part of the general solver by pointint to 
% the dynamics function via a function handle}
% Runs Newtons method to solve a nonilnear system of equations
% Output: Solu matrix of size []
% ***************************************************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *******************************************************v********************************************
function solu = solveNewton(collpts,sysparam,approx)
%% Initialization
IM              = initializeMatrices(collpts);
J_local         = IM.J_local;
J               = IM.J;
ns              = IM.ns;
Elb             = IM.Elb;
Erb             = IM.Erb;
SLb             = IM.SLb;
SRb             = IM.SRb;
R               = IM.R;
LJ_local        = IM.LJ_local;
initial_guess   = sysparam.initial_guess;
BC              = sysparam.BC;
omega           = sysparam.omega;
beta            = sysparam.beta;
N               = collpts.N;
M               = collpts.M;
Neq             = collpts.Neq;
D               = approx.D;
%% Function call
[resid,JacB] = duffingNAE_local(initial_guess,BC,omega,beta,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
solu = initial_guess;
naes = 'NEWTON'; naep = 0.3;
tol = 1e-9; MaxIter = 100; iter = 0; istop = 0;
ti = tic;
%% Newton's Iteration
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
tend = toc(ti);
assignin('base','Newton_time',tend);

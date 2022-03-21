% *******************************************************
% Function solverFsolve() {this function could be a part of the general solver by pointint to 
% the dynamics function via a function handle}
% Runs fsolve to solve a nonilnear system of equations
% Output: Solu matrix of size []
% *******************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *******************************************************
function solu = solveFsolve(collpts,sysparam,approx)
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
BCtype          = sysparam.BCtype;
N               = collpts.N;
M               = collpts.M;
Neq             = collpts.Neq;
D               = approx.D;
%% Fsolve solution
f       = @(X) duffingNAE_local(X,BC,omega,beta,D,N,Neq,ns,M,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype);
options = optimoptions(@fsolve,'Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',true,'StepTolerance',1e-20,'FunctionTolerance',1e-20,'UseParallel',true,'FiniteDifferenceType','central');
% options = optimoptions(@fsolve,'Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',false,'JacobPattern',J_pattern,'StepTolerance',1e-20,'FunctionTolerance',1e-20,'UseParallel',true,'FiniteDifferenceType','central');
ti = tic;
[solu,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);
tend = toc(ti);
assignin('base','Fsolve_time',tend);
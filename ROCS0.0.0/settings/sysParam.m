% *****************************************************************************
% Function sysParam(arg1) [problem specific]
% This function defines the system parameters, boundary conditions, and
% initial guess. 
% Output: struct sysparam
% *****************************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *****************************************************************************
function sysparam = sysParam(collpts)
%% System Parameters
sysparam.epsilon = 1;  %shape parameter
sysparam.omega   = 1;      %rad/sec
sysparam.beta    = 0.9;     %Nonlinearity Coefficient
%% Boundary Conditions 
x10 = 0;        %Initial position
x20 = 0;        %Initial velocity
x1f = 2;        %Final position
x2f = 1;        %Final Velocity  
sysparam.BCtype = 'fixed'; %'fixed','free','P0-Pf','P0-Vf','V0-Pf'
sysparam.BC     = [x10,x20,x1f,x2f]';
%% Initial Guess
N = collpts.N;
M = collpts.M;
% sysparam.initial_guess = [2.35*ones(1*N*M,1); 1.3*ones(1*N*M,1); 2.35*ones(1*N*M,1); 1.3*ones(1*N*M,1)]; 
sysparam.initial_guess = 1*ones( 4*N*M,1);
% sysparam.initial_guess = [1*ones(N*M,1);3.4*ones(N*M,1);1*ones(N*M,1);3.8*ones(N*M,1)];    %works for global collocation
% sysparam.initial_guess = [1*rand(2*N*M,1);.8*rand(2*N*M,1)];
% sysparam.initial_guess = [2.35*ones(2*N*M,1);1.3*ones(2*N*M,1)];
end
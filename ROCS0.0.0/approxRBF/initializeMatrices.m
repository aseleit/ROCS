% ******************************************************************
% Function initializeMatrices(arg1)
% Initializes the matrices used within the solver for speed
% Output: Struct
% ******************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% ******************************************************************
function initialmatrices = initializeMatrices(collpts)
%% Inputs to this function
N               = collpts.N;
M               = collpts.M;
Neq             = collpts.Neq;
%% Initialize Matrices
Jlocal_pattern  = ones(4*N);
J_pattern       = kron(eye(M),Jlocal_pattern);
J_local         = zeros(Neq*N,Neq*N,M);
J               = zeros(Neq*N*M, Neq*N*M); % Local matrices are to be concatenated sequentially
ns              = 0;
Elb             = zeros(1,M);
Erb             = Elb;
SLb             = zeros(Neq,1,M);
SRb             = SLb;
R               = zeros(N,M);
sizeJ_local     = size(J_local);
LJ_local        = sizeJ_local(1);
initialmatrices = struct('Jlocal_pattern', Jlocal_pattern,'J_pattern', J_pattern,'J_local',J_local,...
                         'J',J,'ns',ns,'Elb',Elb,'Erb',Erb, 'SLb',SLb,'SRb',SRb,'R',R,...
                         'sizeJ_local',sizeJ_local,'LJ_local',LJ_local);
end
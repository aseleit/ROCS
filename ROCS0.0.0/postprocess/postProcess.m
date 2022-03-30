% ***************************************************************************************************
% Function postProcess() {this function could be a part of the general solver by pointint to 
% postProcess() extracts the correct values of each state from the input
% solution and order them properly
% Output: struct solu containing ordered propagation of each state
% ***************************************************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% ***************************************************************************************************
function solu = postProcess(xx,collpts)
N       = collpts.N;
M       = collpts.M;
x       = reshape(xx,[4*N,M]);
solu.X1 = x(1:N,:);
solu.X2 = x(N+1:2*N,:);
solu.L1 = x(2*N+1:3*N,:);
solu.L2 = x(3*N+1:4*N,:);
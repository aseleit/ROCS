% ******************************************************
% Function getSolution()
% getSolution() routes the program to the solver chosen
% by the user
% Output: [T,X]
% T: time vector with unrepeated values
% X: struct containing the propagation of states
% *******************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *******************************************************
function [T,X] = getSolution(collpts,sysparam,approx)

xx      = solveFsolve(collpts,sysparam,approx);
solu    = postProcess(xx,collpts);
[T,X]   = removeRepetition(solu,approx.te);



% x = solveNewton(collpts,sysparam,approx);
% compTime_CRBFnewton_end = toc(compTime_CRBFnewton_start);
% SOLU = reshape(solu,[4*N,M]);
% X1_newton = SOLU(1:N,:);
% X2_newton = SOLU(N+1:2*N,:);
% L1_newton = SOLU(2*N+1:3*N,:);
% L2_newton = SOLU(3*N+1:4*N,:);
% 
% X1_newton = X1_newton(:);
% X2_newton = X2_newton(:);
% L1_newton = L1_newton(:);
% L2_newton = L2_newton(:);
% 
% ucrbf_newton = -L2_newton;
% J_newton = 0.5*trapz(ucrbf_newton.^2);




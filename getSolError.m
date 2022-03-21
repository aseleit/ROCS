% ******************************************************
% Function getSolError()
% getSolError() calculates the error between the CRBF
% approximation solution and ode45 solution
% Output: [t,Ex,EL]
% t: time vector
% EX: states error vector of the size [NxM,2]
% EL: costates error vector of the size [NxM,2]
% *******************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *******************************************************
function [t,EX,EL,Xexact] = getSolError(T,X,sysparam)
% IC = [X1(1) X2(2)    3.401297538141612    0.659726264262221];
X1 = X.X1;
X2 = X.X2;
L1 = X.L1;
L2 = X.L2;
IC = [X1(1) X2(1) L1(1) L2(1)];
fexact = @(t,Xexact)duffingDE(Xexact,sysparam.omega,sysparam.beta);
opts = odeset('RelTol',1e-20,'AbsTol',1e-20);
[t, Xexact] = ode23(fexact,T,IC,opts);
X1error = abs(X1' - Xexact(:,1));
% X1error = abs(X1 - x1iclocs);
X2error = abs(X2' - Xexact(:,2));
% X2error = abs(X2 - x2iclocs);
L1error = abs(L1' - Xexact(:,3));
L2error = abs(L2' - Xexact(:,4));
EX = [X1error,X2error];
EL = [L1error, L2error];
% Uerror = abs(L2 - uiclocs);
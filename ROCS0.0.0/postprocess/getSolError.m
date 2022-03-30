% ******************************************************
% Function getSolError()
% getSolError() calculates the error between the CRBF
% approximation solution and ode45 solution
% Output: solutionerror struct
% t: time vector
% X1error: states error vector of the size [NxM,1]
% X2error: costates error vector of the size [NxM,1]
% L1error: 
% *******************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *******************************************************
function solutionerror = getSolError(T,X,sysparam)
% IC = [X1(1) X2(2)    3.401297538141612    0.659726264262221];
X1 = X.X1;
X2 = X.X2;
L1 = X.L1;
L2 = X.L2;
IC = [X1(1) X2(1) L1(1) L2(1)]';
fexact = @(t,Xexact)duffingDE(Xexact,sysparam.omega,sysparam.beta);
opts = odeset('RelTol',1e-14,'AbsTol',1e-14);
[solutionerror.t, solutionerror.Xexact] = ode45(fexact,T,IC,opts);
solutionerror.X1error = abs(X1' - solutionerror.Xexact(:,1));
% X1error = abs(X1 - x1iclocs);
solutionerror.X2error = abs(X2' - solutionerror.Xexact(:,2));
% X2error = abs(X2 - x2iclocs);
solutionerror.L1error = abs(L1' - solutionerror.Xexact(:,3));
solutionerror. L2error = abs(L2' - solutionerror.Xexact(:,4));
% Uerror = abs(L2 - uiclocs);

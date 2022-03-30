format long
clear; close all; clc;
collpts      = collPts();
sysparam     = sysParam(collpts);
approx       = getApprox(collpts,sysparam);
[T,solution] = getSolution(collpts,sysparam,approx);
%% Validation
solutionerror = getSolError(T,solution,sysparam);
%% Plotting
plotSolution(T,solution,collpts,approx,solutionerror)
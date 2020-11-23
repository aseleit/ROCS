clear all;close all;format compact;clc;

[problem,guess]=DuffingOscillator4;          % Fetch the problem definition
options= problem.settings(30);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);

problem.sim.functions=@DuffingOscillator_Dynamics_Sim_Exact;
[ tv1, xv1, uv1 ] = simulateSolution( problem, solution, 'ode23', 0.01 );


figure(1)
plot(tv1,xv1(:,1),'r')
hold on
plot(tv1,xv1(:,2),'b')
plot(tv1,uv1,'k')
axis('tight')
legend('x1','x2','u')



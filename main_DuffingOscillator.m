clear all;close all;format compact;clc;
%new line in this function
%local test one for local collocation
[problem,guess]=DuffingOscillator1;          % Fetch the problem definition
options= problem.settings(30);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
problem.sim.functions=@DuffingOscillator_Dynamics_Sim_Exact;
[ tv1, xv1, uv1 ] = simulateSolution( problem, solution,'ode45');
%% System parameters
epsilon = 10;
tf = 2;
N = 40;%7,9
omega = 1; %rad/sec
x10 = 0;
x20 = 0;
x1f = 5;
x2f = 2;
BCtype = 'fixed';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'
BC = [x10,x20,x1f,x2f]';
%% Approximation matrix
tau = linspace(-1,1,N)';
% tau = lglnodes(N-1);
% tau=flip(tau);
t = mean(tf)/2*tau + mean(tf)/2;
phi = zeros(N,N);
phid = phi;
for i = 1 : N
    phi(i,:)  = rbf0(epsilon,t(i),t);
    phid(i,:) = rbf1(epsilon,t(i),t);
end
D = phid/phi;
%% Solving system of NAE
switch BCtype
    case 'fixed'
        beta = 0.9;
%         initial_guess = [1*ones(2*N,1);3.8*ones(2*N,1)];
%         initial_guess = [10*ones(2*N,1);10*ones(2*N,1)];
        initial_guess = [1*ones(1*N,1);3.5*ones(1*N,1);3.5*ones(1*N,1);3.5*ones(1*N,1)]; 

    case 'free'
        beta = 0.9;
        initial_guess = [1*ones(2*N,1);3.8*ones(2*N,1)]; 
    case 'P0-Pf'
        beta = 0.94;
        initial_guess = [20*ones(2*N,1);1*ones(2*N,1)]; 
    case 'P0-Vf'
        beta = 0.94;
        initial_guess = [20*ones(2*N,1);1*ones(2*N,1)]; 
    case 'V0-Pf'
        beta = 0.97;
        initial_guess = [1*ones(1*N,1);1*ones(1*N,1);0.3*ones(1*N,1);0.3*ones(1*N,1)]; 
end
f = @(X) duffingNAE(X,BC,omega,beta,D,N,BCtype);
options = optimset('Display','iter','Algorithm','Levenberg-Marquardt','Jacobian','off','TolX',1e-14,'TolFun',1e-14);
[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);
X1 = xx(1:N);
X2 = xx(N+1:2*N);
L1 = xx(2*N+1:3*N);
L2 = xx(3*N+1:4*N);
%% Validation
% ticlocs=linspace(solution.T(1,1),solution.tf,N);


x1iclocs = speval(solution,'X',1,t);
x2iclocs = speval(solution,'X',2,t);
uiclocs = speval(solution,'U',1,t);

%ICs from ICLOCS
IC = [x1iclocs(1) x2iclocs(1) -20.3783805876315 1.33559152498820];


fexact = @(t,Xexact)duffingDE(Xexact,omega,beta);
opts = odeset('RelTol',1e-20,'AbsTol',1e-20);
[t, Xexact] = ode45(fexact,t,IC,opts);
X1crbferr = abs(X1 - Xexact(:,1));
X1iclocserr = abs(Xexact(:,1) - x1iclocs);
X1relativeerr = abs(X2 - x2iclocs);

X2crbferr = abs(X2 - Xexact(:,2));
X2iclocserr = abs(Xexact(:,2) - x2iclocs);
% L1error = abs(L1 - Xexact(:,3));
% L2error = abs(L2 - Xexact(:,4));
Uerror = abs(L2 - uiclocs);


%% Plotting
figure(1)
% plot(tv1,xv1(:,1),'r')
plot(t,x1iclocs,'h')

hold on
% plot(tv1,xv1(:,2),'b')
plot(t,x2iclocs,'h')

% plot(tv1,uv1,'k')
plot(t,uiclocs,'h')
axis('tight')
legend('x1','x2','u')

plot(t,X1,'*')
hold on
plot(t,X2,'*')
plot(t,-L2,'*')
axis('tight')
legend('x1','x2','u')

plot(t,Xexact(:,1),'o')
hold on
plot(t,Xexact(:,2),'o')
plot(t,-Xexact(:,4),'o')
axis('tight')
legend('x1','x2','u')
xlabel('time[sec]')
ylabel('Exact states and control')
hold off

figure(3)
plot(t,X1crbferr)
hold on
plot(t,X1iclocserr)
ylabel('X1 abs error')
xlabel('time [sec]')
legend('CRBF error', 'ICLOCS error')

figure(4)
plot(t,X2crbferr)
hold on
plot(t,X2iclocserr)
legend('CRBF error', 'ICLOCS error')
ylabel('X2 abs error')
xlabel('time [sec]')


figure(5)
plot(t,Uerror)


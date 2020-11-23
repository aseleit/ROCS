clear; close all; clc;
%% System parameters
epsilon = 10;
tf = 2;
N = 40;%7,9
omega = 1; %rad/sec
beta = 0.94;
x10 = 0;
x20 = 0;
x1f = 5;
x2f = 2;
BCtype = 'P0-Vf';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'
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
f = @(X) duffingNAE(X,BC,omega,beta,D,N,BCtype);
switch BCtype
    case 'fixed'
        initial_guess = [1*ones(2*N,1);3.8*ones(2*N,1)];
    case 'free'
        initial_guess = [1*ones(2*N,1);3.8*ones(2*N,1)]; 
    case 'P0-Pf'
        initial_guess = [20*ones(2*N,1);1*ones(2*N,1)]; 
    case 'P0-Vf'
        initial_guess = [20*ones(2*N,1);1*ones(2*N,1)]; 
    case 'V0-Pf'
        initial_guess = [1*ones(1*N,1);1*ones(1*N,1);0.3*ones(1*N,1);0.3*ones(1*N,1)]; 
end
options = optimset('Display','iter','Algorithm','levenberg-marquardt','Jacobian','off','TolX',1e-14,'TolFun',1e-14);
[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);
X1 = xx(1:N);
X2 = xx(N+1:2*N);
L1 = xx(2*N+1:3*N);
L2 = xx(3*N:4*N);
%% Validation



%%Plotting
figure(1)
plot(t,X1)
xlabel('time [sec]')
ylabel('x_1')
figure(2)
plot(t,X2)
ylabel('x_2')
xlabel('time [sec]')

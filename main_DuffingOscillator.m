clear all;close all;format compact;clc;
%new line in this function
%local test one for local collocation

%***** Check the final states value in the problem definition *****

BCtype = 'free';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'

% compTime_IClOCS_start = tic;
% [problem,guess]=DuffingOscillator2;          % Fetch the problem definition
% options= problem.settings(40);                  % Get options and solver settings
% [solution,MRHistory]=solveMyProblem( problem,guess,options);
% problem.sim.functions=@DuffingOscillator_Dynamics_Sim_Exact;
% [ tv1, xv1, uv1 ] = simulateSolution( problem, solution,'ode45');
% compTime_IClOCS_end = toc(compTime_IClOCS_start);
% EPSILON = [1 10 100 1000 1e4 1e6];
EPSILON = 100;
for count = 1 : length(EPSILON)
%% System parameters
epsilon = EPSILON(count);
% epsilon = 1;
tic
tf = 2;
N = 40;% Number of nodes
% epsilon = (N-1)/4/tf;
omega = 1; %rad/sec
x10 = 0;
x20 = 0;
x1f = 2;
x2f = 1;
BC = [x10,x20,x1f,x2f]';
%% Approximation matrix
tau = linspace(-1,1,N)';
% tau = lglnodes(N-1);
% tau = flip(tau);
t = mean(tf)/2*tau + mean(tf)/2;
phi = zeros(N,N);
phid = phi;
qt = phi;
qdt = phi;
for i = 1 : N
    phi(i,:)  = rbf0(epsilon,t(i),t);
    phid(i,:) = rbf1(epsilon,t(i),t);
    qt(i,:) = phi(i,:) / sum(phi(i,:));
    qdt(i,:) = (phid(i,:)*sum(phi(i,:)) - phi(i,:)*sum(phid(i,:)))/(sum(phi(i,:)))^2;
end
D = phid/phi;
% D = qdt/qt;

%% Solving system of NAE
switch BCtype
    case 'fixed'
        beta = .9;
        %         initial_guess = [3.8*ones(2*N,1);1*ones(2*N,1)];
%                   initial_guess = zeros(4*N,1);
        initial_guess = [2.35*ones(2*N,1);1.3*ones(2*N,1)];
%                 initial_guess = [1.78*ones(N,1);1*ones(N,1);0*ones(N,1);0*ones(N,1)];
        %         initial_guess = [0*ones(1*N,1);10*ones(1*N,1);0*ones(1*N,1);0*ones(1*N,1)];
        
    case 'free'
        beta = 0.9;
        %         initial_guess = [1*ones(2*N,1);3.8*ones(2*N,1)];
        %         initial_guess = [2.35*ones(2*N,1);1.3*ones(2*N,1)];
        initial_guess = [10*ones(1*N,1);10*ones(1*N,1);0*ones(1*N,1);0*ones(1*N,1)];
        
        
    case 'P0-Pf'
        beta = 0.94;
        % initial_guess = [20*ones(2*N,1);1*ones(2*N,1)];
        % works for Jacobian Off 
%          initial_guess = [10*ones(1*N,1);10*ones(1*N,1);0*ones(1*N,1);0*ones(1*N,1)]; 
        initial_guess = [2.3*ones(2*N,1);2.3*ones(2*N,1)];

        
    case 'P0-Vf'
        beta = 0.94;
        initial_guess = [20*ones(2*N,1);1*ones(2*N,1)];
    case 'V0-Pf'
        beta = 0.97;
        initial_guess = [1*ones(1*N,1);1*ones(1*N,1);0.3*ones(1*N,1);0.3*ones(1*N,1)];
end
f = @(X) duffingNAE(X,BC,omega,beta,D,N,BCtype);
options = optimset('Display','iter','Algorithm','Levenberg-Marquardt','Jacobian','on','TolX',1e-20,'TolFun',1e-20);
compTime_CRBFfsolve_start = tic;
[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);
compTime_CRBFfsolve_end = toc(compTime_CRBFfsolve_start); % computational time elapsed
X1 = xx(1:N);
X2 = xx(N+1:2*N);
L1 = xx(2*N+1:3*N);
L2 = xx(3*N+1:4*N);
U = -L2;
J = 0.5*trapz(U.^2);
JJ = 0.5*(U'*U);
%% Shooting method
disp('********** SHOOTING METHOD *************')
shooting_guess = [-0.4 -0.4]' ;
compTime_shooting_start = tic;
f2 = @(XXX) duffingShooting(XXX,omega,beta,BC,t,BCtype);
% Find the initial Costates
options2 = optimset('Display','iter','Algorithm','Levenberg-Marquardt','TolX',1e-20,'TolFun',1e-20);
Li_shooting = fsolve(f2,shooting_guess,options2); 
% Propagate with the initial costates
fshooting = @(t,x_shooting)duffingDE(x_shooting,omega,beta);
options_shooting = odeset('RelTol',1e-6,'AbsTol',1e-6);
IC_shooting = [x10; x20; Li_shooting];
[t, x_shooting] = ode45(fshooting,t,IC_shooting,options_shooting); 
compTime_shooting_end = toc(compTime_shooting_start);
%[x_shooting,fval2,exitflag2,output2]
u_shooting = -x_shooting(:,4);
J_shooting =  0.5*trapz(u_shooting.^2);
%% Newton's method

switch BCtype 
    case 'fixed'
        initial_guess = [2.35*ones(2*N,1);1.3*ones(2*N,1)];
        [resid,JacB] = duffingNAE(initial_guess_newton,BC,omega,beta,D,N,BCtype);
        solu = initial_guess_newton;
    case 'free'
        initial_guess_newton = [0*ones(2*N,1); 0*ones(2*N,1)];
        [resid,JacB] = duffingNAE(initial_guess_newton,BC,omega,beta,D,N,BCtype);
        solu = initial_guess_newton;
    case 'P0-Pf'
        
    case 'P0-Vf'
        
    case 'V0-Pf'
        
end

naes = 'NEWTON'; naep = 0.3;
tol = 1e-9; MaxIter = 100; iter = 0; istop = 0;
compTime_CRBFnewton_start = tic;
while(istop==0)
    if(iter>MaxIter)
        disp('Max iteration reached for the nonliner solver')
        istop = 1;
    end
    dsolu = nae_su(resid,JacB,solu,initial_guess,naes,naep,iter); % compute incremental solution   
    solu = solu + dsolu; 
    [resid,JacB] = duffingNAE(solu,BC,omega,beta,D,N,BCtype);
    % determine convergence
    if(norm(resid)<tol)
        istop = 1;
    end
    iter = iter + 1;
    fprintf('Iteration = %d\n',iter)
end
compTime_CRBFnewton_end = toc(compTime_CRBFnewton_start);
X1_newton = solu(1:N);
X2_newton = solu(N+1:2*N);
L1_newton = solu(2*N+1:3*N);
L2_newton = solu(3*N+1:4*N);
ucrbf_newton = -L2_newton;
J_newton = 0.5*trapz(ucrbf_newton.^2);
%% Validation
ticlocs=linspace(solution.T(1,1),solution.tf,N);
x1iclocs = speval(solution,'X',1,t);
x2iclocs = speval(solution,'X',2,t);
uiclocs = speval(solution,'U',1,t);
%jiclocs = 0.5*trapz(uiclocs.^2);
jjiclocs = 0.5*(uiclocs'*uiclocs);
%ICs from ICLOCS
% IC = [x10 x20 -0 0];
IC = [X1(1), X2(1), L1(1), L2(1)];

fexact = @(t,Xexact)duffingDE(Xexact,omega,beta);
opts = odeset('RelTol',1e-20,'AbsTol',1e-20);
tic
[t, Xexact] = ode45(fexact,t,IC,opts);
toc
X1exact = Xexact(:,1);
X2exact = Xexact(:,2);
L1exact = Xexact(:,3);
L2exact = Xexact(:,4);

X1crbferr = abs(X1 - X1exact);
% X1iclocserr = abs(Xexact(:,1) - x1iclocs);
% X1relativeerr = abs(X2 - x2iclocs);

X2crbferr = abs(X2 - X2exact);
% X2iclocserr = abs(Xexact(:,2) - x2iclocs);
% L1error = abs(L1 - Xexact(:,3));
L2error = abs(L2 - L2exact);
% Uerror = abs(U) - abs(uiclocs);

x1RMSE(count,1) =  sqrt(mean(X1crbferr).^2/N);
x2RMSE(count,1) = sqrt(mean(X2crbferr).^2/N);
uRMSE(count,1) = sqrt(mean(L2error).^2/N);
end
%% Plotting
switch BCtype
    case "fixed"
        figure(1)
        plot(t,X1,'b*')
        hold on
        plot(t,X2,'b+')
        plot(t,-L2,'bx')        
        plot(t,Xexact(:,1),'ro')
        plot(t,Xexact(:,2),'rd')
        plot(t,-Xexact(:,4),'rs')
        legend({'x1(CRBF)','x2(CRBF)','u(CRBF)','x1(exact)','x2(exact)','u(exact)'},'Location','northwest')
        axis('tight')
        xlabel('Time [sec]')
        ylabel('States and control')
        hold off
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(1), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\states.eps','epsc2');
%         
%         
        figure(2)
        plot(t,X1crbferr,'b')
        hold on
        ylabel('Absolute error in x_1 ')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(2), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\x1err.eps','epsc2');
%         
        
        figure(3)
        plot(t,X2crbferr,'b')
        hold on
        ylabel('Absolute error in x_2')
        xlabel('Time [sec]')
        grid on
        set(gcf,'renderer','painters') 
%         saveas(figure(3), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\x2err.eps','epsc2');
%         
        figure(4)
        plot(t,L2error,'b')
        ylabel('Absolute error in u')
        xlabel('Time [sec]')
        grid on
        set(gcf,'renderer','painters') 
%         saveas(figure(4), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\L2err.eps','epsc2');
%         
        figure(5)
        semilogy(t,X1crbferr,'-or')
        hold on
        semilogy(t,X2crbferr,'-db')
        semilogy(t,L2error,'-sk')
        ylim([-3e-3 3e-3])
        ylabel('Absolute error in states and control')
        xlabel('Time [sec]')
        legend('x_1','x_2','u')
        grid on
        set(gcf,'renderer','painters') 
%         saveas(figure(5), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\allerr.eps','epsc2');
%    
%         figure(6)
%         bar(categorical({'J_{CRBF}','J_{ICLOCS2}'}),[J,jiclocs])
%         ylabel('Energy Cost')
%         xlabel('Approximation Method')

        figure(6)
        plot(t,U)
        hold on
        plot(t,uiclocs)
        plot(t,-L2exact)
        legend('CRBF','ICLOCS','Exact')
    case "free"
        figure(1)
        plot(t,X1,'b*')
        hold on
        plot(t,X2,'b+')
        plot(t,-L2,'bx')        
        plot(t,Xexact(:,1),'ro')
        plot(t,Xexact(:,2),'rd')
        plot(t,-Xexact(:,4),'rs')
        legend({'x1(CRBF)','x2(CRBF)','u(CRBF)','x1(exact)','x2(exact)','u(exact)'},'Location','northwest')
        axis('tight')
        xlabel('Time [sec]')
        ylabel('States and control')
        hold off
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(1), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\states.eps','epsc2');
%         
        
      figure(2)
        plot(t,X1crbferr,'b')
        hold on
        ylabel('Absolute error in x_1 ')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(2), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\x1err.eps','epsc2');
%         
        
       figure(3)
        plot(t,X2crbferr,'b')
        hold on
        ylabel('Absolute error in x_2')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(3), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\x2err.eps','epsc2');
%         
       figure(4)
        plot(t,L2error,'b')
        ylabel('Absolute error in u')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(4), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\L2err.eps','epsc2');
%     
        figure(5)
        plot(t,X1crbferr,'-or')
        hold on
        plot(t,X2crbferr,'-db')
        plot(t,L2error,'-sk')
        ylabel('Absolute error in states and control')
        xlabel('Time [sec]')
        legend('x_1','x_2','u')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(5), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\allerr.eps','epsc2');
%            figure(6)
        plot(t,U)
        hold on
        plot(t,uiclocs)
        plot(t,-L2exact)
        legend('CRBF','ICLOCS','Exact')

    case "P0-Pf"
       figure(1)
        plot(t,X1,'b*')
        hold on
        plot(t,X2,'b+')
        plot(t,-L2,'bx')         
        plot(t,Xexact(:,1),'ro')
        plot(t,Xexact(:,2),'rd')
        plot(t,-Xexact(:,4),'rs')
        legend({'x1(CRBF)','x2(CRBF)','u(CRBF)','x1(exact)','x2(exact)','u(exact)'},'Location','northwest')
        axis('tight')
        xlabel('Time [sec]')
        ylabel('States and control')
        hold off
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(1), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\states.eps','epsc2');
%         
        
      figure(2)
        plot(t,X1crbferr,'b')
        hold on
        ylabel('Absolute error in x_1 ')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(2), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\x1err.eps','epsc2');
%         
        
        figure(3)
        plot(t,X2crbferr,'b')
        hold on
        ylabel('Absolute error in x_2')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(3), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\x2err.eps','epsc2');
%         
        figure(4)
        plot(t,L2error,'b')
        ylabel('Absolute error in u')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(4), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\L2err.eps','epsc2');
%         
        figure(5)
        plot(t,X1crbferr,'-or')
        hold on
        plot(t,X2crbferr,'-db')
        plot(t,L2error,'-sk')
        ylabel('Absolute error in states and control')
        xlabel('Time [sec]')
        legend('x_1','x_2','u')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(5), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\allerr.eps','epsc2');
%    
        figure(6)
        plot(t,U)
        hold on
        plot(t,uiclocs)
        plot(t,-L2exact)
        plot(t,u_shooting)
        legend('CRBF','ICLOCS','Exact','Shooting')
        
        
        
       case "P0-Vf"
        figure(1)
        plot(t,X1,'b*')
        hold on
        plot(t,X2,'b+')
        plot(t,-L2,'bx')     
        plot(t,Xexact(:,1),'ro')
        plot(t,Xexact(:,2),'rd')
        plot(t,-Xexact(:,4),'rs')
        legend({'x1(CRBF)','x2(CRBF)','u(CRBF)','x1(exact)','x2(exact)','u(exact)'},'Location','northwest')
        axis('tight')
        xlabel('Time [sec]')
        ylabel('States and control')
        hold off
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(1), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\states.eps','epsc2');
%         
%         
        figure(2)
        plot(t,X1crbferr,'b')
        hold on
        ylabel('Absolute error in x_1 ')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(2), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\x1err.eps','epsc2');
%         
        
        figure(3)
        plot(t,X2crbferr,'b')
        hold on
        ylabel('Absolute error in x_2')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(3), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\x2err.eps','epsc2');
%         
        figure(4)
        plot(t,L2error,'b')
        ylabel('Absolute error in u')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(4), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\L2err.eps','epsc2');
%         
        figure(5)
        semilogy(t,X1crbferr,'-or')
        hold on
        semilogy(t,X2crbferr,'-db')
        semilogy(t,L2error,'-sk')
        axis([0 2 -3e-3 3e-3])
        ylabel('Absolute error in states and control')
        xlabel('Time [sec]')
        legend('x_1','x_2','u')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(5), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\allerr.eps','epsc2');
%         
        case "V0-Pf"
        figure(1)
        plot(t,X1,'b*')
        hold on
        plot(t,X2,'b+')
        plot(t,-L2,'bx')     
        plot(t,Xexact(:,1),'ro')
        plot(t,Xexact(:,2),'rd')
        plot(t,-Xexact(:,4),'rs')
        legend({'x1(CRBF)','x2(CRBF)','u(CRBF)','x1(exact)','x2(exact)','u(exact)'},'Location','northwest')
        axis('tight')
        xlabel('Time [sec]')
        ylabel('States and control')
        hold off
        grid on
        set(gcf,'renderer','painters') 
        saveas(figure(1), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\states.eps','epsc2');
        
        
        figure(2)
        plot(t,X1crbferr,'b')
        hold on
        ylabel('Absolute error in x_1 ')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(2), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\x1err.eps','epsc2');
%         
%         
        figure(3)
        plot(t,X2crbferr,'b')
        hold on
        ylabel('Absolute error in x_2')
        xlabel('Time [sec]')
        grid on
        set(gcf,'renderer','painters') 
        saveas(figure(3), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\x2err.eps','epsc2');
        
        figure(4)
        plot(t,L2error,'b')
        ylabel('Absolute error in u')
        xlabel('Time [sec]')
        grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(4), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\L2err.eps','epsc2');
%         
        figure(5)
        semilogy(t,X1crbferr,'-or')
        axis([0 2 -3e-3 3e-3])
        hold on
        semilogy(t,X2crbferr,'-db')
        semilogy(t,L2error,'-sk')
        ylabel('Absolute error in states and control')
        xlabel('Time [sec]')
        legend('x_1','x_2','u')
        grid on
        set(gcf,'renderer','painters') 
%         saveas(figure(5), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\allerr.eps','epsc2');
end
%% Plotting RMSE
% switch BCtype
%     case 'fixed'
%         figure(10)
%         semilogx(EPSILON,x1RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_1');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(10), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\x1RMSE.eps','epsc2');
%         
%         figure(11)
%         semilogx(EPSILON,x2RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_2');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(11), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\x2RMSE.eps','epsc2');
%         
%         figure(12)
%         semilogx(EPSILON,uRMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in u');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(12), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\uRMSE.eps','epsc2');
%         
%         figure(13)
%         semilogx(EPSILON,x1RMSE,'-r')
%         hold on
%         semilogx(EPSILON,x2RMSE,'--b')
%         semilogx(EPSILON,uRMSE,'-.k')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in states and control');
%         legend('x_1','x_2','u')
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(13), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_fixed\allRMSE.eps','epsc2');
%         
%         
%     case 'free'
%         figure(10)
%         semilogx(EPSILON,x1RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_1');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(10), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\x1RMSE.eps','epsc2');
%         
%         figure(11)
%         semilogx(EPSILON,x2RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_2');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(11), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\x2RMSE.eps','epsc2');
%         
%         figure(12)
%         semilogx(EPSILON,uRMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in u');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(12), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\uRMSE.eps','epsc2');
%         
%         figure(13)
%         semilogx(EPSILON,x1RMSE,'-r')
%         hold on
%         semilogx(EPSILON,x2RMSE,'--b')
%         semilogx(EPSILON,uRMSE,'-.k')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in states and control');
%         legend('x_1','x_2','u')
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(13), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_free\allRMSE.eps','epsc2');
%         
%     case 'P0-Pf'
%         figure(10)
%         semilogx(EPSILON,x1RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_1');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(10), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\x1RMSE.eps','epsc2');
%         
%         figure(11)
%         semilogx(EPSILON,x2RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_2');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(11), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\x2RMSE.eps','epsc2');
%         
%         figure(12)
%         semilogx(EPSILON,uRMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in u');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(12), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\uRMSE.eps','epsc2');
%         
%         
%         figure(13)
%         semilogx(EPSILON,x1RMSE,'-r')
%         hold on
%         semilogx(EPSILON,x2RMSE,'--b')
%         semilogx(EPSILON,uRMSE,'-.k')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in states and control');
%         legend('x_1','x_2','u')
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(13), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PP\allRMSE.eps','epsc2');
%         
%     case 'P0-Vf'
%         figure(10)
%         semilogx(EPSILON,x1RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_1');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(10), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\x1RMSE.eps','epsc2');
%         
%         figure(11)
%         semilogx(EPSILON,x2RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_2');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(11), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\x2RMSE.eps','epsc2');
%         
%         figure(12)
%         semilogx(EPSILON,uRMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in u');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(12), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\uRMSE.eps','epsc2');
%         
%         
%         figure(13)
%         semilogx(EPSILON,x1RMSE,'-r')
%         hold on
%         semilogx(EPSILON,x2RMSE,'--b')
%         semilogx(EPSILON,uRMSE,'-.k')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in states and control');
%         legend('x_1','x_2','u')
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(13), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_PV\allRMSE.eps','epsc2');
%         
%     case 'V0-Pf'
%         figure(10)
%         semilogx(EPSILON,x1RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_1');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(10), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\x1RMSE.eps','epsc2');
%         
%         figure(11)
%         semilogx(EPSILON,x2RMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in x_2');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(11), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\x2RMSE.eps','epsc2');
%         
%         figure(12)
%         semilogx(EPSILON,uRMSE,'b')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in u');
%         grid on
%         set(gcf,'renderer','painters') 
%         saveas(figure(12), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\uRMSE.eps','epsc2');
%         
%         
%         figure(13)
%         semilogx(EPSILON,x1RMSE,'-r')
%         hold on
%         semilogx(EPSILON,x2RMSE,'--b')
%         semilogx(EPSILON,uRMSE,'-.k')
%         xlabel('Shape parameter c');
%         ylabel('RMSE in states and control');
%         legend('x_1','x_2','u')
%         grid onh˙
%         set(gcf,'renderer','painters') 
%         saveas(figure(13), 'D:\OneDrive - Knights - University of Central Florida\ARMY RESEARCH LAB\CRBF OCP Paper\Nonlinear Dynamics Journal Template\figures\duffing_VP\allRMSE.eps','epsc2');
% end



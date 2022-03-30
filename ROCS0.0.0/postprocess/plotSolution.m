function plotSolution(T,solution,collpts,approx,err)
X1 = solution.X1;
X2 = solution.X2;
L1 = solution.L1;
L2 = solution.L2;
te = approx.te;
t = te(:);
% figure(1)
plot(t,X1,'*-')
hold on
plot(t,X2,'+-')
% plot(t,-L2,'*')
axis('tight')
% legend('x1','x2','u')
grid on
plot(t,err.Xexact(:,1),'o')
hold on
plot(t,err.Xexact(:,2),'o')
% plot(t,-Xexact(:,4),'o')
axis('tight')
legend('x1','x2','u','x1_exact','x2_exact','u_exact')
xlabel('time[sec]')
ylabel('Exact states and control')
for pp = 1:collpts.M-1
    tline = te(:,end,pp);
    xline(tline,'linewidth',2);
end
hold off
figure(2)
plot(T,L1)
hold on
plot(T,L2)
grid on
legend('lambda_1','lambda_2')
for pp = 1:collpts.M-1
    tline = te(:,end,pp);
    xline(tline,'linewidth',2);
end
hold off

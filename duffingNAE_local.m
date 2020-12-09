function R = duffingNAE_local(X,BC,omega,beta,D,N,M)
%% State, Costates and BCs 
% this needs review x1 = X(1:N)
x1 = X(1:N*M);
x1 = reshape(x1,[N,M]);
x2 = X(N*M+1 : 2*N*M);
x2 = reshape(x2,[N,M]);
L1 = X(2*N*M+1 : 3*N*M);
L1 = reshape(L1,[N,M]);
L2 = X(3*N*M+1 : 4*N*M);
L2 = reshape(L2,[N,M]);

% for k = 1:M
%     Elb(k) = (k-1)*LA+1;    %left boundary of each segment
%     Erb(k) = Elb(k)+LA-1;   %right boundary of each segment
%     for ii = 1 : 4
%         Slb(ii,k) = Elb(k)+ (ii-1)*n; %left boundary of each state ii for each element k
%         Srb(ii,k) = Elb(k) + ii*n-1;
%     end
% end



% x1 = [X(1:N); X(4*N+1:5*N)];
% x2 = [X(N+1:2*N); X(5*N+1:6*N)];
% L1 = [X(2*N+1:3*N); X(6*N+1:7*N)];
% L2 = [X(3*N+1:4*N); X(7*N+1:8*N)];
% x1 = reshape(x1,[N,M]);
% x2 = reshape(x2,[N,M]);
% L1 = reshape(L1,[N,M]);
% L2 = reshape(L2,[N,M]);

x10 = BC(1);
x20 = BC(2);
x1f = BC(3);
x2f = BC(4);

R = zeros(N,M);
for k = 1 : M
    R(1:N,k) = D(:,:,k)*x1(:,k) - x2(:,k);% Nx1 = NxN * Nx1 - Nx1;
    R(N+1:2*N,k) = D(:,:,k)*x2(:,k) + omega^2*x1(:,k) + beta*x1(:,k).^3 + L2(:,k);
    R(2*N+1:3*N,k) = D(:,:,k)*L1(:,k) - L2(:,k).*(omega^2 + 3*beta*x1(:,k).^2);
    R(3*N+1:4*N,k) = D(:,:,k)*L2(:,k) + L1(:,k);
end
% Continuity conditions
% if M > 1
    for kk = 2 : M
%         R(1,kk) = -R(N,kk-1);
%         R(N+1,kk) = -R(2*N,kk-1);
%         R(2*N+1,kk) = -R(3*N,kk-1);
%         R(3*N+1,kk) = -R(4*N,kk-1);

        R(1,kk) = x1(1,kk) - x1(N,kk-1);
        R(N+1,kk) = x2(1,kk) - x2(end,kk-1);
        R(2*N+1,kk) = L1(1,kk) - L1(end,kk-1);
        R(3*N+1,kk) = L2(1,kk) - L2(end,kk-1);

    end
% end
% Boundary conditions
R(N+1,1) = x1(1,1) - x10; % IC 1 in element 1
R(2*N+1,1) = x2(1,1) - x20; % IC 2 in element 1
R(2*N,end) = x1(end,end) - x1f; % FC 1 in last element
R(3*N,end) = x2(end,end) - x2f; % FC 2 in last element
% Function output
R = R(:);  %[all state segment1; all states segment2; ... ; all states last segment];
%[x2 all sements; x2 all segments; L1 all segments; L2 all segments]
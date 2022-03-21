% *****************************************************************
% Function getApprox(arg1,arg2)
% This function generates the approximation matrices/tensors
% getApprox() generates a tensor if the number of segments M > 1
% D:      Differential operator of the size [NxNxM]
% phi:    RBF gram matrix of the size [NxNxM]
% phid:   RBF gram matrix derivative of the size [NxNxM]
% te:     Tensor of segmented time
% Output: struct approx with three fields
% ******************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% ******************************************************************
function approx = getApprox(collpts,sysparam)
N       = collpts.N;
M       = collpts.M;
tf      = collpts.tf;
epsilon = sysparam.epsilon;
D       = zeros(N,N,M);
phi     = zeros(N,N,M);
phid    = zeros(N,N,M);
for k = 1 : M
    t0e            = (k-1)*tf/M; % Initial time of each segment (element)
    tfe            = k*tf/M;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N);
    %     tinterp(:,:,k) = linspace(t0e,tfe,5*n);
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
        %         phi_interp(i,:,k) = rbf0(c,te(:,i,k),tinterp(:,:,k));
    end
    D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end
approx.D    = D;
approx.phid = phid;
approx.phi  = phi;
approx.te = te;
end

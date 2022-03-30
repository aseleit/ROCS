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
epsilon                  = sysparam.epsilon;
approx.te                = getTimeSegments(collpts);
approx.r                 = getRadiiMatrix(collpts,approx.te);
[approx.phi,approx.phid] = getRBF(approx.r,epsilon,collpts.RBFtype);
approx.D                 = approx.phid/approx.phi;
end

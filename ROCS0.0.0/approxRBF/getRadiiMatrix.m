% *************************************************************************
% Function getRadiiMatrix()
% getRadiiMatrix() generates the radii matrix to be plugged in an RBF
% The distance from each point in the time vector is subtracted from the 
% whole time vector. The process is repeated to generate every row of the 
% r matrix. The process is also repeated to generate a matrix for every
% time segment of index k.
% Output: r tensor of the size NxNxM
% *************************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *************************************************************************
function r = getRadiiMatrix(collpts,te)
N       = collpts.N;
M       = collpts.M;
r = zeros(N,N,M);
for k = 1 : M
    for i = 1 : N
            r(i,:,k) = (te(:,i,k) - te(:,:,k));
    end
end
% *************************************************************************
% Function getTimeSegments()
% getTimeSegments() divides the time domain into M segments 
% Each segment contains N collocation points
% Output: te 3D tensor with a series of collocation points in each segment 
% te(:,:,index of segment)
% *************************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *************************************************************************
function te = getTimeSegments(collpts)
N       = collpts.N;
M       = collpts.M;
tf      = collpts.tf;
te      = zeros(1,N,M);          % Matrix Initializtion for speed
for k = 1 : M
    t0e            = (k-1)*tf/M; % Initial time of each segment (element)
    tfe            = k*tf/M;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N);
end

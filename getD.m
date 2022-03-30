function D = getD(phid,phi,collpts)
D = zeros(collpts.N,collpts.N,collpts.M);
for k = 1 : collpts.M 
    D(:,:,k) = phid(:,:,k)/phi(:,:,k);
end
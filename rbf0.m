function phi = rbf0(c,t,ti)

r = ti - t;
R = abs(r);

h = 1/c*R;

F = sqrt(h.^2+1);

phi = F + R.^5;


%CIMQ
% phi = 1./F + R.^5;
% phi = sqrt(1 + (c*R).^2); %multiquadric



end
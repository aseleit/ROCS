function phid = rbf1(c,ti,t)


r = ti - t;
% r = t - ti;
R = abs(r);

h = 1/c*R;

% Coupled Multiquadric
phid = (h./c./sqrt(h.^2+1)+ 5*R.^4).*sign(r);
% phid = (h./c./sqrt(h.^2+1)-2*c^2*R.*exp(-c^2*R.^2)).*sign(r);

% Coupled Gaussian
% phid = (-2*h./c.*exp(-h.^2) + 5*R.^4).*sign(r);
% CIMQ
% phid = (h./c./(h.^2+1).^(3/2) + 5*R.^4).*sign(r);
% Multiquadric
% phid = c^2*R./(sqrt(1 + (c*R).^2)).*sign(r);

% Gaussian
% phid = -2*c^2*R.*exp(-c^2*R.^2).*sign(r);
end
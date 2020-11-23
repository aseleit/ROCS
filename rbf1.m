function phid = rbf1(c,t,ti)


r = t - ti;
R = abs(r);

h = 1/c*R;

phid = (h./c./sqrt(h.^2+1) + 5*R.^4).*sign(r);

% CIMQ
% phid = (h./c./(h.^2+1).^(3/2) + 5*R.^4).*sign(r);
% Multiquadric
% phid = c^2*R./(sqrt(1 + (c*R).^2));

end
function phi = rbf0(c,ti,t)

r = ti - t;
% r = t - ti;
R = abs(r);

h = 1/c*R;

F = sqrt(h.^2+1);

% Coupled Multiquadric
phi = F + R.^5 ;
% phi = F + exp(-c^2*R.^2);




% Coupled Gaussian
% phi = exp(-h.^2) + R.^5 ;
%CIMQ
% phi = 1./F + R.^5;
% phi = sqrt(1 + (c*R).^2); %multiquadric

%gaussian
% phi = exp(-c^2*R.^2);


%%NRBF
%phi(i,:) = qt(i,:) / sum(qt(i,:));
%phid(i,:) = (qdt(i,:)*sum(qt(i,:)) - qt(i,:)*sum(qdt(i,:)))/(sum(qt(i,:)))^2;

end
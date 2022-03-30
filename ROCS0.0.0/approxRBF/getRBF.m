% *************************************************************************
% Function getRBF()
% getRBF() generates the RBF gram matrix phi, its derivative phid and its
% integral phiI.
% Output: [phi,phid,phiI]
% *************************************************************************
% Ahmed Seleit, 2022, Aerospace Engineering, UCF
% *************************************************************************
function [phi,phid] = getRBF(R,c,RBFtype)
h = R/c;
F = sqrt(h.^2+1);
switch RBFtype
    case 'MCQ'  % Coupled Multiquadric
        phi  = F + R.^5 ;
        phid = (h./c./sqrt(h.^2+1)+ 5*R.^4).*sign(R); %review sign(r)
    case 'CGA'  % Coupled Gaussian
        phi  = exp(-h.^2) + R.^5 ;
        phid = (-2*h./c.*exp(-h.^2) + 5*R.^4).*sign(R);
    case 'CIMQ' % Coupled Inverse Multiquadric
        phi  = 1./F + R.^5;
        phid = (h./c./(h.^2+1).^(3/2) + 5*R.^4).*sign(R);
    case 'R1'
        phi = h;
    case 'R3'
        phi = h.^3;
    case 'TPS2'
        I = (h > 0);
        phi(I) = h(I).^2.*log(h(I));
    case 'Q'
        phi = 1 + h.^2;
    case 'MQ'
        phi = sqrt(1 + h.^2);
    case 'IMQ'
        phi = 1./sqrt(1 + h.^2);
    case 'IQ'
        phi = 1./(1 + h.^2);
    case 'GS'
        phi = exp(-h.^2);
    case 'CP_C0'
        I = (h < 1);
        phi(I) = (1 - h(I)).^2;
    case 'CP_C2'
        I = (h < 1);
        phi(I) = (1 - h(I)).^4.*(4*h(I) + 1);
    case 'CP_C4'
        I = (h < 1);
        phi(I) = (1 - h(I)).^6.*(35/3*h(I).^2 + 6*h(I) + 1);
    case 'CP_C6'
        I = (h < 1);
        phi(I) = (1 - h(I)).^8.*(32*h(I).^3 + 25*h(I).^2 + 8*h(I) + 1);
    case 'CTPS_C0'
        I = (h < 1);
        phi(I) = (1 - h(I)).^5;
    case 'CTPS_C1'
        I = (h < 1 & h > 0);
        phi(I) = 1 + 80/3*h(I).^2 - 40*h(I).^3 + 15*h(I).^4 - 8/3*h(I).^5 + 20*h(I).^2.*log(h(I));
        phi(h == 0) = 1;
    case 'CTPS_C2a'
        I = (h < 1 & h > 0);
        phi(I) = 1 - 30*h(I).^2 - 10*h(I).^3 + 45*h(I).^4 - 6*h(I).^5 - 60*h(I).^3.*log(h(I));
        phi(h == 0) = 1;
    case 'CTPS_C2b'
        I = (h < 1 & h > 0);
        phi(I) = 1 - 20*h(I).^2 + 80*h(I).^3 - 45*h(I).^4 -16*h(I).^5 + 60*h(I).^4.*log(h(I));
        phi(h == 0) = 1;
end
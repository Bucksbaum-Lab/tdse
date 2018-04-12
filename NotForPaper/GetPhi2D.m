function [phi1, phi2, Emin, deltaEgrid] = GetPhi2D(W11, W22, W12, psi10, psi20, dx, mm, nn, xx)
hbar = 0.6582;

Kmax = hbar^2/mm/2*(pi/dx)^2/2;

Emin = min(min([min(W11), min(W22), min(W12)]));
deltaEgrid = max(max([max(W11), max(W22), max(W12)]))+Kmax-Emin;

Hnorm1 = 2*(W11 - (deltaEgrid/2+Emin))/deltaEgrid;
Hnorm2 = 2*(W22 - (deltaEgrid/2+Emin))/deltaEgrid;
Hnorm12 = 2*W12/deltaEgrid;

KEfactor = -2*4*hbar^2/2/mm/deltaEgrid;

phi1 = zeros(size(psi10,1),size(psi10,2),nn);
phi2 = zeros(size(psi10,1),size(psi10,2),nn);

phi1(:,:,1) = psi10;
phi1(:,:,2) = -1i*(KEfactor*del2(psi10,xx,xx)+Hnorm1.*psi10+Hnorm12.*psi20);

%sum(sum(conj(phi1(:,:,2)).*psi10))

phi2(:,:,1) = psi20;
phi2(:,:,2) = -1i*(KEfactor*del2(psi20,xx,xx)+Hnorm2.*psi20+Hnorm12.*psi10);

for ii = 3:nn
    
    phi1(:,:,ii) = -2*1i*(KEfactor*del2(phi1(:,:,ii-1),xx,xx)+Hnorm1.*phi1(:,:,ii-1)+Hnorm12.*phi2(:,:,ii-1))+phi1(:,:,ii-2);
    phi2(:,:,ii) = -2*1i*(KEfactor*del2(phi2(:,:,ii-1),xx,xx)+Hnorm2.*phi2(:,:,ii-1)+Hnorm12.*phi1(:,:,ii-1))+phi2(:,:,ii-2);
    
end
%sum(sum(conj(phi1(:,:,3)).*psi10))
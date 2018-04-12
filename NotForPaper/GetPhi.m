function [phi1, phi2, Emin, deltaEgrid] = GetPhi(V1, V2, V12, psi10, psi20, dx, mm, nn, xx)
hbar = 0.6582;

Kmax = hbar^2/mm/2*(pi/dx)^2/4;

Emin = min([min(V1), min(V2), min(V12)]);
deltaEgrid = max([max(V1), max(V2), max(V12)])+Kmax-Emin;

Hnorm1 = 2*(V1 - (deltaEgrid/2+Emin))/deltaEgrid;
Hnorm2 = 2*(V2 - (deltaEgrid/2+Emin))/deltaEgrid;
Hnorm12 = 2*V12/deltaEgrid;

KEfactor = -8*hbar^2/2/mm/deltaEgrid;

phi1 = zeros(nn,length(psi10));
phi2 = zeros(nn,length(psi10));

phi1(1,:) = psi10;
phi1(2,:) = -1i*(KEfactor*del2(psi10,xx)+Hnorm1.*psi10+Hnorm12.*psi20);

phi2(1,:) = psi20;
phi2(2,:) = -1i*(KEfactor*del2(psi20,xx)+Hnorm2.*psi20+Hnorm12.*psi10);

for ii = 3:nn
    
    phi1(ii,:) = -2*1i*(KEfactor*del2(phi1(ii-1,:),xx)+Hnorm1.*phi1(ii-1,:)+Hnorm12.*phi2(ii-1,:))+phi1(ii-2,:);
    phi2(ii,:) = -2*1i*(KEfactor*del2(phi2(ii-1,:),xx)+Hnorm2.*phi2(ii-1,:)+Hnorm12.*phi1(ii-1,:))+phi2(ii-2,:);
    
end
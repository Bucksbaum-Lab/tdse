function [phi1, phi2, Emin, deltaEgrid] = OldGetPhi(V1, V2, V12, psi10, psi20, dx, NN, mm, nn)
hbar = 0.6582;

kk = pi*linspace(-1/dx, 1/dx, NN);
kk = fftshift(kk);

T = hbar^2/mm/2*kk.^2/4;


Emin = min([min(V1), min(V2), min(V12)]);
deltaEgrid = max([max(V1), max(V2), max(V12)])+max(T)-Emin;

T = 2*T/deltaEgrid;

Hnorm1 = 2*(V1 - (deltaEgrid/2+Emin))/deltaEgrid;
Hnorm2 = 2*(V2 - (deltaEgrid/2+Emin))/deltaEgrid;
Hnorm12 = 2*V12/deltaEgrid;

phi1 = zeros(nn,length(psi10));
phi2 = zeros(nn,length(psi10));

phi1(1,:) = psi10;
phi1(2,:) = -1i*(ifft(T.*fft(psi10))+Hnorm1.*psi10+Hnorm12.*psi20);

phi2(1,:) = psi20;
phi2(2,:) = -1i*(ifft(T.*fft(psi20))+Hnorm2.*psi20+Hnorm12.*psi10);

for ii = 3:nn
    
    phi1(ii,:) = -2*1i*(ifft(T.*fft(phi1(ii-1,:)))+Hnorm1.*phi1(ii-1,:)+Hnorm12.*phi2(ii-1,:))+phi1(ii-2,:);
    phi2(ii,:) = -2*1i*(ifft(T.*fft(phi2(ii-1,:)))+Hnorm2.*phi2(ii-1,:)+Hnorm12.*phi1(ii-1,:))+phi2(ii-2,:);
    
end
%Figuring out the momentum operator

close all
clear all
clc

LL = 50;
NN = 2^10;
xx = linspace(-LL, LL, NN);
dx = LL/NN;
mm = 1;
hbar = 1;
a=1;

kk = linspace(-1/dx, 1/dx, NN);
kk = fftshift(kk);
KEUU = -pi^2/4*(kk.^2);

wavefunction = sin(pi*xx/a);

momentum = ifft(KEUU.*fft(wavefunction));

figure()
plot(xx, wavefunction, 'LineWidth', 1)
hold on
plot(xx, -(pi/a)^2*sin(pi*xx/a), 'LineWidth', 3)
plot(xx, real(momentum), 'LineWidth', .5)
plot(xx, imag(momentum), 'LineWidth', 1)
axis([-10,10,-20,20])
legend('wave', 'deriviative', 'KE', 'imag')


[x1,x2] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

wavefunction = sin(pi*x1/a).*sin(pi*x2/a);
KEUU = -pi^2/4*(k1.^2+k2.^2);
momentum = ifft2(KEUU.*fft2(wavefunction));

derivative = -2*(pi/a)^2*sin(pi*x1/a).*sin(pi*x2/a);

figure()
surf(x1,x2, wavefunction)
colorbar
axis([-1,1,-1,1])
caxis([-20,20])
%hold on
figure()
subplot(2,1,1)
surf(x1,x2, derivative)
colorbar
axis([-1,1,-1,1])
caxis([-20,20])
subplot(2,1,2)
surf(x1, x2, real(momentum))
colorbar
axis([-1,1,-1,1])
caxis([-20,20])
%figure()
%contour(x1, x2, imag(momentum),4)
%colorbar
%axis([-1,1,-1,1])
%caxis([-20,20])
%axis([-10,10,-20,20])
%legend('wave', 'deriviative', 'KE', 'imag')

%}
clear all
close all
clc

%% Discritize the wavefunction

LL = 100;
NN = 2^15;
xx = linspace(-LL, LL-2*LL/NN, NN); %Angstrom
dt = 1/100; %femtoseconds
totalTime = 10; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

%% Constants

kappa = 1; %eV/Angstrom
lambda = 0; %eV/Angstrom
mm = 10*0.0096; %amu
hbar = 0.6582;
lighteV = 10; %eV
nu = 10; %eV/Angstrom

omegaLight = lighteV/hbar;

%% Hamiltonians

%W11 = kappa/2*xx.*xx;
W11 = kappa*xx.^2/2;
W12 = lambda*xx;
W22 = kappa*xx.^2+10;

%% Initial wavefunctions

Dwave1 = exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%+2*sqrt(sqrt(kappa*mm)/hbar)*xx.*exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%Dwave1 = sin(xx);
Dwave1 = Dwave1/sqrt(sum(sum(Dwave1.*Dwave1)));
Dwave2 = zeros(size(Dwave1));
%Dwave2 = exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%Dwave2 = Dwave2/sqrt(sum(sum(Dwave2.*Dwave2)));
%Dwave1 = zeros(size(Dwave2));

'del2'
sum((-4*del2(Dwave1,xx)*hbar^2/2/mm+W11.*Dwave1).*conj(Dwave1))
plot(xx,-4*del2(Dwave1,xx)*hbar^2/2/mm)

'FFT'
kk = pi*linspace(-1/dx, 1/dx, NN);
kk = fftshift(kk);
T = hbar^2/mm/2*kk.^2/4;
sum((ifft(T.*fft(Dwave1))+W11.*Dwave1).*conj(Dwave1))

plot(xx,-4*del2(Dwave1,xx)*hbar^2/2/mm, xx, ifft(T.*fft(Dwave1)))
legend('del2','FFT')

'actual'
(hbar*sqrt(kappa/mm))/2
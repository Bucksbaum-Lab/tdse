%check time

%% Discritize the wavefunction

LL = 20;
NN = 2^9;
xx = linspace(-LL, LL-2*LL/NN, NN);
dx = LL/NN;
mm = 25;
dt = 1/1000;
tNN = 1;

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[XX,YY] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

%% Constants

kappa = 1;
lambda = 0;
x0 = 2.5;
omegay = 2;
my = 5;
hbar = 1;

%% Hamiltonians

W11 = kappa*XX;
W12 = lambda*YY;
W22 = -kappa*XX;

[Vp, Vm, WD, A, Ainv] = makeW(W11, W12, W22);
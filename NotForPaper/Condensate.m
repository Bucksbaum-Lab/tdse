%Split Step 2D

clear all
close all

%Step one - Discritize the wavefunction

LL = 10;
NN = 2^10;
xx = linspace(-LL, LL, NN);
dx = LL/NN;
mm = 1;
dt = -1i/1000;
tNN = 2000;
hbar = 1;
omega = 1;

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[x1,x2] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

%wavefunction = exp(-mm*omega/2/hbar*(x1.^2/100+x2.^2/100));
%wavefunction = x1.*x2.*exp(-mm*omega/(2*hbar)*(x1.^2+x2.^2));
wavefunction = rand(size(x1));
wavefunction = wavefunction/sqrt(sum(sum(wavefunction.*wavefunction)));

KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2));

VV = mm*omega^2/2*(x1.^2+x2.^2);
VB = -exp(-sqrt(x1.^2+x2.^2)/50);
%VB = 0;

VU = exp(-(1i*VV+VB)/hbar*dt);

figure()
surf(x1, x2, abs(wavefunction))
shading flat
drawnow
%figure()
%surf(x1, x2, angle(wavefunction))
%shading flat
%drawnow

%Step two - propagate

for ll = 1:tNN

    tic
    
    if mod(ll,200) == 0
        figure()
        surf(x1, x2, abs(wavefunction))
        shading flat
        drawnow
        %bool = abs(wavefunction) > 0.00001;
        %Angle = angle(wavefunction);
        %Angle = bool.*Angle;
        %figure()
        %surf(x1, x2, Angle)
        %shading flat
        %drawnow
    end
    
    wavefunction = VU.*ifft2(KEUU.*fft2(wavefunction));
    wavefunction = wavefunction/sqrt(sum(sum(real(wavefunction.*conj(wavefunction)))));
    
    if mod(ll,20) == 0
        (tNN-ll)*toc
        'seconds left'
    end
end

        figure()
        surf(x1, x2, abs(wavefunction))
        shading flat
        drawnow
        %bool = abs(wavefunction) > 0.00001;
        %Angle = angle(wavefunction);
        %Angle = bool.*Angle;
        %figure()
        %surf(x1, x2, Angle)
        %shading flat
        %drawnow
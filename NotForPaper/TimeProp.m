close all
clear all

%numbers

%angstroms
L = 50;
fs = L/2^9;
m = 1;
N = L/fs;
hbar = 1.54e-34*6.24e18*1e15;
%fs
ts = 1/2^7;

%x space, momentum space, and time
x = -L/2:fs:L/2-fs;
k = fftshift((-N/2:N/2-1)/L);

%build hamiltonian
KE = exp(-hbar^2/(2*m)*(-(2*pi*k).^2/(1i*hbar)*ts));
PE = exp(-x/(1i*hbar)*ts);

%build initial wavefunction
psi0 = exp(-x.^2/2);
wavefunction = psi0/sqrt(sum(psi0.^2*2*L*fs));

nmax = 200;

Mag(nmax+1) = struct('cdata',[],'colormap',[]);
Angle(nmax+1) = struct('cdata',[],'colormap',[]);

figure()
plot(x,wavefunction.*conj(wavefunction))
title('propagated')
drawnow
Mag(1) = getframe(gcf);

figure()
plot(x,angle(wavefunction))
title('phase')
drawnow
Angle(1) = getframe(gcf);

for nn = 1:nmax
    
    close all
    
    momentum = fft(wavefunction);
    wavefunction = PE.*(ifft(KE.*momentum));
    
    figure()
    plot(x,wavefunction.*conj(wavefunction))
    title('propagated')
    drawnow
    Mag(nn+1) = getframe(gcf);
    
    figure()
    plot(x,unwrap(angle(wavefunction)))
    title('phase')
    drawnow
    Angle(nn+1) = getframe(gcf);
    
    wavefunction = wavefunction./sqrt(real(sum(wavefunction.*conj(wavefunction)*2*L*fs)));
    
end

%wave = psi0*exp(-1i*.5*t(nn))+psi1*exp(-1i*1.5*(t(nn)))+psi2*exp(-1i*2.5*(t(nn)));
%wave = wave/sqrt(sum(wave.*conj(wave)*2*L*fs));

%{
figure()
subplot(2,2,1)
plot(x,real(wavefunction))
title('propagated, real')
subplot(2,2,2)
plot(x,real(wave))
title('exact, real')

subplot(2,2,3)
plot(x,imag(wavefunction))
title('propagated, imag')
subplot(2,2,4)
plot(x,imag(wave))
title('exact, imag')
%}
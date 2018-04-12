function Energy = FourierMethod(wavefunction, Potential, k, m, fs, omega)

hbar = 1.54e-34;
%hbar=1;

momentum = fftshift(fft(wavefunction));
Kinetic = -(2*pi*k).^2.*momentum;
Der2 = ifft(ifftshift(Kinetic));

Energy = (Potential.*wavefunction-hbar^2/(2*m)*Der2);


secondDer = (wavefunction-2*circshift(wavefunction,[0,1])+circshift(wavefunction,[0,2]))/(fs*sqrt(hbar/m/omega))^2;


Energyder = Potential.*wavefunction-hbar^2/(2*m)*secondDer;


%Energy(30000)/wavefunction(30000)/hbar
%Energyder(30000)/wavefunction(30000)/hbar
%Der2(30000)/secondDer(30000)

%isreal(hbar^2/(2*m)*Der2.*conj(Der2))

%Der2(2e3)/secondDer(2e3)

%{
figure()
subplot(2,2,3)
plot(real(hbar^2/(2*m)*Der2.*conj(Der2)))
title('Fourier Method')
subplot(2,2,4)
plot(real(hbar^2/(2*m)*secondDer.*conj(secondDer)))
title('derivative')
subplot(2,2,1)
plot(real(Energy.*conj(Energy))/hbar^2)
title('energy')
subplot(2,2,2)
plot(k,real(momentum.*conj(momentum)))
title('momentum')
%}
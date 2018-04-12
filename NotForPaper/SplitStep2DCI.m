%Split Step 2D

clear all
close all
clc

%Step one - Discritize the wavefunction

LL = 100;
NN = 2^10;
xx = linspace(-LL, LL, NN);
dx = LL/NN;
mm = 25;
dt = 1/10;
tNN = 100000;
hbar = 1;
omega = 1;

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[x1,x2] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

wavefunctionE = exp(-mm/2/hbar*((x1-50).^2/100+(x2-50).^2/100));
wavefunctionE = wavefunctionE/sqrt(sum(sum(wavefunctionE.*wavefunctionE)));
wavefunctionG = zeros(size(wavefunctionE));

zlimit = max(max(abs(wavefunctionE))).*1.1;

wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);

f = figure();
plotThatStuff(wavefunctionE, wavefunctionG, x1, x2, zlimit, LL)
drawnow
wavemovie(1) = getframe(f);
close all

KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2));

%Step two - set up the CI

delta = 1;
lambda = 1;

W11 = -delta*x1;
W22 = delta*x1;
W12 = lambda*x2;

VG = -(.5*(W11+W22)+sqrt((.5*(-W22+W11)).^2+W12.^2));
VE = -(.5*(W11+W22)-sqrt((.5*(-W22+W11)).^2+W12.^2));

UG = exp(-1i*VG/hbar*dt);
UE = exp(-1i*VE/hbar*dt);

%Adiabatic to diabatic transform

AA11 = cos(.5*atan(2*W12./(W11-W22)));
AA12 = sin(.5*atan(2*W12./(W11-W22)));
AA11(x1==0&x2==0) = cos(.5*atan(lambda/(delta)));
AA12(x1==0&x2==0) = sin(.5*atan(lambda/(delta)));

DwaveG = AA11.*wavefunctionG - AA12.*wavefunctionE;
DwaveE = AA12.*wavefunctionG + AA11.*wavefunctionE;

%Step three - propagate

time = 0;

tic

for ll = 1:tNN
 
    if mod(ll,25) == 0
        
        f = figure();
        plotThatStuff(wavefunctionE, wavefunctionG, x1, x2, zlimit, LL)
        
        wavemovie(ll/25+1) = getframe(f);
        close all

        time = toc + time;
        tic
        averagetime = time/ll;
        
        disp([num2str((tNN-ll)*averagetime) ' seconds left'])
    end
    
    DwaveE = ifft2(KEUU.*fft2(DwaveE));
    DwaveG = UG.*ifft2(KEUU.*fft2(DwaveG));
        
    wavefunctionG = AA11.*DwaveG + AA12.*DwaveE;
    wavefunctionE = -AA12.*DwaveG + AA11.*DwaveE;
    
    wavefunctionG = UG.*wavefunctionG;
    wavefunctionE = UE.*wavefunctionE;
    
    wavefunctionG = wavefunctionG/sqrt(sum(sum(real(wavefunctionG.*conj(wavefunctionG)))));
    wavefunctionE = wavefunctionE/sqrt(sum(sum(real(wavefunctionE.*conj(wavefunctionE)))));
    
    DwaveG = AA11.*wavefunctionG - AA12.*wavefunctionE;
    DwaveE = AA12.*wavefunctionG + AA11.*wavefunctionE;
    
end

x=linspace(0,2001,200);
y=sin(x);
sound(y)

t = datetime;

save(['C:\Chelsea\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])
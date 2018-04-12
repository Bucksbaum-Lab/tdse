%Split Step 2D

clear all
close all
clc

%Step one - Discritize the wavefunction

m = 0.58;
A = 3.0;
D = 5.0;
sigma = 0.3;
sigma1 = 0.75;
omega0 = 39.14e13;
omega1 = 7.83e13;



kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[x1,x2] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

wavefunctionE = exp(-mm/2/hbar*((x1+25).^2/1+(x2+25).^2/1));
wavefunctionE = wavefunctionE/sqrt(sum(sum(wavefunctionE.*wavefunctionE)));
wavefunctionG = zeros(size(wavefunctionE));

zlimit = 0.04;

wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);

plotThatStuff(wavefunctionE, wavefunctionG, x1, x2, zlimit, LL)

drawnow

wavemovie(1) = getframe(gcf);

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

for ll = 1:tNN

    tic
    
    if mod(ll,25) == 0
        %plotThatStuff(wavefunctionE, wavefunctionG, x1, x2, zlimit, LL)
        %drawnow
        f = figure()
        
        subplot(2,2,1)
        surf(x1, x2, abs(wavefunctionE))
        shading flat
        axis([-LL,LL,-LL,LL,0,zlimit])
        %caxis([0,zlimit])
        xlabel('x1')
        ylabel('x2')

        bool = abs(wavefunctionE) > 0.00001;
        Angle = angle(wavefunctionE);
        Angle = bool.*Angle;

        subplot(2,2,2)
        surf(x1, x2, Angle)
        shading flat
        axis([-LL,LL,-LL,LL,0,zlimit])
        %caxis([0,zlimit])
        xlabel('x1')
        ylabel('x2')

        subplot(2,2,3)
        surf(x1, x2, abs(wavefunctionG))
        shading flat
        axis([-LL,LL,-LL,LL,0,zlimit])
        %caxis([0,zlimit])
        xlabel('x1')
        ylabel('x2')

        bool = abs(wavefunctionG) > 0.00001;
        Angle = angle(wavefunctionG);
        Angle = bool.*Angle;

        subplot(2,2,4)
        surf(x1, x2, Angle)
        shading flat
        axis([-LL,LL,-LL,LL,0,zlimit])
        %caxis([0,zlimit])
        xlabel('x1')
        ylabel('x2')

        drawnow
        wavemovie(ll/25+1) = f;
        close all
    end
    
    DwaveE = ifft2(KEUU.*fft2(DwaveE));
    DwaveG = UG.*ifft2(KEUU.*fft2(DwaveG));
        
    wavefunctionG = AA11.*DwaveG + AA12.*DwaveE;
    wavefunctionE = -AA12.*DwaveG + AA11.*DwaveE;
    
    wavefunctionG = UG.*wavefunctionG;
    wavefunctionE = UE.*wavefunctionE;
    
    wavefunctionG = wavefunctionG/sqrt(sum(sum(real(wavefunctionG.*conj(wavefunctionG)))));
    wavefunctionE = wavefunctionE/sqrt(sum(sum(real(wavefunctionE.*conj(wavefunctionE)))));
        
    if mod(ll,5) == 0
        (tNN-ll)*toc
        'seconds left'
    end
end

x=linspace(0,2001,200);
y=sin(x);
sound(y)

%save('C:\Chelsea\TDSE\OutputOct6-2015.mat')
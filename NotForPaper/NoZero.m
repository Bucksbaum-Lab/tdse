%Split Step 2D

%{
_     _
|Vp  0|
|0  Vm|
-     -
%}

clear all
close all
clc
%delete(gcp('nocreate'))
%parpool(8)

%% Discritize the wavefunction

LL = 40;
NN = 2^8;
xx = linspace(-LL, LL, NN); %Angstrom
dt = 1/6000; %femtoseconds
totalTime = 2.2; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[XX,YY] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

%% Constants

kappa = 2; %eV/Angstrom
lambda = .5; %eV/Angstrom
mm = 28*0.0096; %amu
hbar = 0.658;
lighteV = 5; %eV
nu = 0; %eV/Angstrom
xposition = 5;
width = 1;
timeposition = 1;
timewidth = 5;

omegaLight = lighteV/hbar;

%% Hamiltonians

%W11 = kappa*XX.^2/100;
%W11 = kappa/2*(XX.^2+YY.^2);
W22 = -kappa*XX;
%W11 = -kappa*(XX-50);
%.*(1-tanh(XX-50))+kappa*YY.^2;
W12 = lambda*YY;
%W12 = lambda*exp(-XX.^2).*YY - nu*sin(omegaLight*0);
%W22 = kappa/2*(XX.^2+YY.^2)+1;
W11 = kappa*XX;
%W22 = kappa*(50+XX);
%.*(1-tanh(-50-XX))+kappa*YY.^2;

%KEUU = exp(-1i*hbar*pi^2/2/mm*dt*(k1.^2+k2.^2)/8);
%KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2)/2);
KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2)/2);

C = (W11-W22).^2;
A = exp(-1i/hbar*(W11+W22)*dt/2);
B = (W22-W11);

D = sqrt(4*W12.*conj(W12)+C);

U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
U12 = A.*1i.*sin(D*dt/2/hbar)./D*(-2).*W12;

%U11(D == 0) = A(D==0).*(1+1i*dt/2/hbar.*B(D==0));
%U22(D == 0) = A(D==0).*(1-1i*dt/2/hbar.*B(D==0));
%U12(D == 0) = A(D==0).*1i*dt/2/hbar*(-2).*W12(D==0);

%% Initial wavefunctions

%Dwave1 = exp(-mm/2/hbar*((XX+2).^2+(YY).^2)/.1);
Dwave2 = exp(-(((XX+xposition)/width).^2+(YY/width).^2));
Dwave2 = Dwave2/sqrt(sum(sum(Dwave2.*Dwave2)));
Dwave1 = zeros(size(Dwave2));

D0 = Dwave2;

zlimit = max(max(abs(Dwave2)))^2.*1.1;
zlimit2 = max(sum(Dwave2.*conj(Dwave2),1))*1.1;

%% Stuff for figures

wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
wavemovie2 = wavemovie;
anglemovie = wavemovie;
wavemovieA = wavemovie;
wavemovie2A = wavemovie;
anglemovieA = wavemovie;
jj = 1000;

%% Propagate

time = 0;
conservation = zeros(tNN,1);
Excited = conservation;
Ground = conservation;
electricField = conservation;

disp('Coming up with time estimate...')

tic

%% Adiabatic to diabatic
AA = zeros(NN,NN,2,2);
Ainv = AA;
WD = AA;

for nn = 1:NN
    for mm = 1:NN
        
        [AA(nn,mm,:,:), WD(nn,mm,:,:), Ainv(nn,mm,:,:)] = eig([W11(nn,mm), W12(nn,mm); W12(nn,mm), W22(nn, mm)]);
                
    end
end

DiabaticToAdiabatic11 = AA(:,:,1,1);
DiabaticToAdiabatic12 = AA(:,:,1,2);
DiabaticToAdiabatic22 = AA(:,:,2,2);
DiabaticToAdiabatic21 = AA(:,:,2,1);

Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic21.*Dwave2;
Awave2 = DiabaticToAdiabatic22.*Dwave2 + DiabaticToAdiabatic12.*Dwave1;

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,jj) == 1
                
        f = figure();
        subplot(2,1,1)
        surf(XX, YY, Dwave2.*conj(Dwave2))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['D excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX, YY, Dwave1.*conj(Dwave1))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['D ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovie((ll-1)/jj+1) = getframe(f);
        
        g = figure();
        subplot(2,2,3)
        plot(XX(1,:), sum(Dwave1.*conj(Dwave1),1))
        axis([-LL LL 0 zlimit2])
        title(['D ground X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,4)
        plot(YY(:,1), sum(Dwave1.*conj(Dwave1),2))
        axis([-LL LL 0 zlimit2])
        title(['D ground Y sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,1)
        plot(XX(1,:), sum(Dwave2.*conj(Dwave2),1))
        axis([-LL LL 0 zlimit2])
        title(['D excited X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,2)
        plot(YY(:,1), sum(Dwave2.*conj(Dwave2),2))
        axis([-LL LL 0 zlimit2])
        title(['D excited Y sums ' num2str(ll*dt) ' fs'])
        wavemovie2((ll-1)/jj+1) = getframe(g);
        
        bool = abs(Dwave2) > 0.0001;
        Angle2 = angle(Dwave2);
        Angle2 = bool.*Angle2;
        
        bool = abs(Dwave1) > 0.0001;
        Angle1 = angle(Dwave1);
        Angle1 = bool.*Angle1;
        
        h = figure();
        subplot(2,1,1)
        pcolor(XX, YY, Angle2.*conj(Angle2))
        axis([-LL LL -LL LL])
        caxis([0,2*pi])
        xlabel('x')
        ylabel('y')
        title(['D excited ' num2str(ll*dt) ' fs'])
        shading flat
        colorbar
        subplot(2,1,2)
        pcolor(XX, YY, Angle1.*conj(Angle1))
        axis([-LL LL -LL LL])
        caxis([0,2*pi])
        xlabel('x')
        ylabel('y')
        title(['D ground ' num2str(ll*dt) ' fs'])
        shading flat
        colorbar
        anglemovie((ll-1)/jj+1) = getframe(h);
        
        F = figure();
        subplot(2,1,1)
        surf(XX, YY, Awave2.*conj(Awave2))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['A excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX, YY, Awave1.*conj(Awave1))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['A ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovieA((ll-1)/jj+1) = getframe(F);
        
        G = figure();
        subplot(2,2,3)
        plot(XX(1,:), sum(Awave1.*conj(Awave1),1))
        axis([-LL LL 0 zlimit2])
        title(['A ground X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,4)
        plot(YY(:,1), sum(Awave1.*conj(Awave1),2))
        axis([-LL LL 0 zlimit2])
        title(['A ground Y sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,1)
        plot(XX(1,:), sum(Awave2.*conj(Awave2),1))
        axis([-LL LL 0 zlimit2])
        title(['A excited X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,2)
        plot(YY(:,1), sum(Awave2.*conj(Awave2),2))
        axis([-LL LL 0 zlimit2])
        title(['A excited Y sums ' num2str(ll*dt) ' fs'])
        wavemovie2A((ll-1)/jj+1) = getframe(G);
        
        bool = abs(Awave2) > 0.0001;
        Angle2 = angle(Awave2);
        Angle2 = bool.*Angle2;
        
        bool = abs(Awave1) > 0.0001;
        Angle1 = angle(Awave1);
        Angle1 = bool.*Angle1;
        
        H = figure();
        subplot(2,1,1)
        pcolor(XX, YY, Angle2.*conj(Angle2))
        axis([-LL LL -LL LL])
        caxis([0,2*pi])
        xlabel('x')
        ylabel('y')
        title(['A excited ' num2str(ll*dt) ' fs'])
        shading flat
        colorbar
        subplot(2,1,2)
        pcolor(XX, YY, Angle1.*conj(Angle1))
        axis([-LL LL -LL LL])
        caxis([0,2*pi])
        xlabel('x')
        ylabel('y')
        title(['A ground ' num2str(ll*dt) ' fs'])
        shading flat
        colorbar
        anglemovieA((ll-1)/jj+1) = getframe(H);
        
        close all
        
        if ll ~= 1
            seconds = (tNN-ll)/jj*toc;
            tic
            hours = floor(seconds/60/60);
            minutes = floor((seconds - hours*60*60)/60);
            seconds = seconds - hours*60*60 - minutes*60;
            disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])
        else
            second = toc;
        end

    end
    
    %% Hamiltonians
    
    if dt*ll > timeposition && dt*ll < timeposition+timewidth
        electricField(ll) = nu*sin(omegaLight*(dt*ll-timeposition));
        W12 = lambda*YY+electricField(ll);
    else
        electricField(ll) = 0;
        W12 = lambda*YY;
    end
        
    D = sqrt(4*W12.*conj(W12)+C);

    U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
    U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
    U12 = A.*1i.*sin(D*dt/2/hbar)./D*(-2).*W12;

    %U11(D == 0) = A(D==0).*(1+1i*dt/2/hbar.*B(D==0));
    %U22(D == 0) = A(D==0).*(1-1i*dt/2/hbar.*B(D==0));
    %U12(D == 0) = A(D==0).*1i*dt/2/hbar*(-2).*W12(D==0);

    %% Propagation
    
    Dwave1 = ifft2(KEUU.*fft2(Dwave1));
    Dwave2 = ifft2(KEUU.*fft2(Dwave2));
    
    Dwave1 = U11.*Dwave1 + U12.*Dwave2;
    Dwave2 = U22.*Dwave2 + U12.*Dwave1;
    
    Dwave1 = ifft2(KEUU.*fft2(Dwave1));
    Dwave2 = ifft2(KEUU.*fft2(Dwave2));

    %Dwave1 = Dwave1.*(.5-(XX/LL/2).^2-(YY/LL/2).^2)*2;
    %exp(-(XX/50).^2-(YY/50).^2);
    %Dwave2 = Dwave2.*(.5-(XX/LL/2).^2-(YY/LL/2).^2)*2;
    %exp(-(XX/50).^2-(YY/50).^2);
    
    conservation(ll) = sum(sum(Dwave2.*conj(Dwave2)))+sum(sum(Dwave1.*conj(Dwave1)));
    
    %Dwave1 = Dwave1.*exp(-(XX/100).^2-(YY/100).^2);
    %Dwave2 = Dwave2.*exp(-(XX/100).^2-(YY/100).^2);
    
    %Dwave1 = Dwave1/conservation(ll);
    %Dwave2 = Dwave2/conservation(ll);
    
    Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic21.*Dwave2;
    Awave2 = DiabaticToAdiabatic22.*Dwave2 + DiabaticToAdiabatic12.*Dwave1;
    
    Excited(ll) = sum(sum(Awave2.*conj(Awave2)));
    Ground(ll) = sum(sum(Awave1.*conj(Awave1)));
    
    if ll == 1
    
        seconds = (tNN-ll)*(toc+(1-jj)/jj*second);
        hours = floor(seconds/60/60);
        minutes = floor((seconds - hours*60*60)/60);
        seconds = seconds - hours*60*60 - minutes*60;
        disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])
        
    end
end
%}
%% Done :)

t = datetime;

%save('C:\Chelsea\TDSE_Output\A0eV.mat')
%save('C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\I-k2-l_5-n1-eV5-m28-sin.mat', '-v7.3')
%save(['C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])
%save(['C:\Chelsea\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])

x=linspace(0,2001,200);
y=sin(x);
sound(y)
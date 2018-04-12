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

LL = 50;
NN = 2^9;
xx = linspace(-LL, LL-2*LL/NN, NN); %Angstrom
dt = 1/5000; %femtoseconds
totalTime = 1; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[XX,YY] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

%% Constants

kappa = 1; %eV/Angstrom
lambda = 0.001; %eV/Angstrom
mm = 1*0.0096; %amu
hbar = 0.658;
lighteV = 4; %eV
nu = 0.01; %eV/Angstrom
xposition = 0.5;
width = 1;

omegaLight = lighteV/hbar;

%% Hamiltonians

%W11 = kappa*XX.^2/100;
%W11 = kappa/2*(XX.^2);
%.*(1-tanh(XX-50))
W11 = kappa*(XX)+10*YY.^2+10*XX.^2;
%W11 = -kappa*(5-XX).*(1-tanh(-5+XX))+500+kappa*YY.^2;
W12 = lambda*YY;
%W22 = kappa/2*(XX.^2+YY.^2)+1;
W22 = -kappa*(XX)+10*YY.^2+10*XX.^2;
%W22 = -kappa*(XX-5).*(1-tanh(XX-5))+kappa*YY.^2;

%KEUU = exp(-1i*hbar*pi^2/2/mm*dt*(k1.^2+k2.^2)/8);
%KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2)/2);
KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2));

C = (W11-W22).^2;
A = exp(-1i/hbar*(W11+W22)*dt/2);
B = (W22-W11);

D = sqrt(4*W12.*conj(W12)+C);

U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
U12 = A.*1i.*sin(D*dt/2/hbar)./D*(-2).*W12;

U11(D == 0) = A(D==0).*(1+1i*dt/2/hbar.*B(D==0));
U22(D == 0) = A(D==0).*(1-1i*dt/2/hbar.*B(D==0));
U12(D == 0) = A(D==0).*1i*dt/2/hbar*(-2).*W12(D==0);

DiabaticToAdiabatic11 = cos(0.5*atan(W12./0.5*(W11-W22)));
DiabaticToAdiabatic12 = sin(0.5*atan(W12./0.5*(W11-W22)));

Dipole = - nu*sin(omegaLight*0)+zeros(size(U11));

UDipole11 = cos(2*abs(Dipole)*dt/2/hbar);
UDipole12 = -1i*sin(2*abs(Dipole)*dt/2/hbar).*sign(Dipole);

%% Initial wavefunctions

%Dwave1 = exp(-mm/2/hbar*((XX+2).^2+(YY).^2)/.1);
Awave2 = exp(-(((XX+xposition)/width).^2+(YY/width).^2));
Awave2 = Awave2/sqrt(sum(sum(Awave2.*Awave2)));
Awave1 = zeros(size(Awave2));

%D0 = Dwave2;

Dwave1 = DiabaticToAdiabatic11.*Awave1 - DiabaticToAdiabatic12.*Awave2;
Dwave2 = DiabaticToAdiabatic12.*Awave1 + DiabaticToAdiabatic11.*Awave2;

zlimit = max(max(abs(Dwave2)))^2.*1.1;
zlimit2 = max(sum(Dwave2.*conj(Dwave2),1))*1.1;

%% Stuff for figures

wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
wavemovie2(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
jj = 100;

%% Propagate

time = 0;
conservation = zeros(tNN,1);
Excited = conservation;
Ground = conservation;

disp('Coming up with time estimate...')

tic

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,jj) == 1
        
        %{
        f = figure();
        subplot(2,1,1)
        surf(XX, YY, Dwave2.*conj(Dwave2))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX, YY, Dwave1.*conj(Dwave1))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovie((ll-1)/jj+1) = getframe(f);
        
        g = figure();
        subplot(2,2,1)
        plot(XX(1,:), sum(Dwave1.*conj(Dwave1),1))
        axis([-LL LL 0 zlimit2])
        title(['ground X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,2)
        plot(YY(:,1), sum(Dwave1.*conj(Dwave1),2))
        axis([-LL LL 0 zlimit2])
        title(['ground Y sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,3)
        plot(XX(1,:), sum(Dwave2.*conj(Dwave2),1))
        axis([-LL LL 0 zlimit2])
        title(['excited X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,4)
        plot(YY(:,1), sum(Dwave2.*conj(Dwave2),2))
        axis([-LL LL 0 zlimit2])
        title(['excited Y sums ' num2str(ll*dt) ' fs'])
        wavemovie2((ll-1)/jj+1) = getframe(g);
        close all
        %}
        
        f = figure();
        subplot(2,1,1)
        surf(XX, YY, Awave2.*conj(Awave2))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX, YY, Awave1.*conj(Awave1))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovie((ll-1)/jj+1) = getframe(f);
        
        g = figure();
        subplot(2,2,1)
        plot(XX(1,:), sum(Awave1.*conj(Awave1),1))
        axis([-LL LL 0 zlimit2])
        title(['ground X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,2)
        plot(YY(:,1), sum(Awave1.*conj(Awave1),2))
        axis([-LL LL 0 zlimit2])
        title(['ground Y sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,3)
        plot(XX(1,:), sum(Awave2.*conj(Awave2),1))
        axis([-LL LL 0 zlimit2])
        title(['excited X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,4)
        plot(YY(:,1), sum(Awave2.*conj(Awave2),2))
        axis([-LL LL 0 zlimit2])
        title(['excited Y sums ' num2str(ll*dt) ' fs'])
        wavemovie2((ll-1)/jj+1) = getframe(g);
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
    
    Dipole =  nu*sin(omegaLight*(ll-1)*dt)+zeros(size(Dwave1));

    UDipole11 = cos(2*abs(Dipole)*dt/2/hbar);
    UDipole12 = -1i*sin(2*abs(Dipole)*dt/2/hbar).*sign(Dipole);

    %% Propagation
    
    Dwave1 = ifft2(KEUU.*fft2(Dwave1));
    Dwave2 = ifft2(KEUU.*fft2(Dwave2));
    
    Dwave1 = U11.*Dwave1 + U12.*Dwave2;
    Dwave2 = U22.*Dwave2 + U12.*Dwave1;
    
    Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic12.*Dwave2;
    Awave2 = -DiabaticToAdiabatic12.*Dwave1 + DiabaticToAdiabatic11.*Dwave2;
    
    Awave1 = UDipole11.*Awave1 + UDipole12.*Awave2;
    Awave2 = UDipole11.*Awave2 + UDipole12.*Awave1;
    
    Dwave1 = DiabaticToAdiabatic11.*Awave1 - DiabaticToAdiabatic12.*Awave2;
    Dwave2 = DiabaticToAdiabatic12.*Awave1 + DiabaticToAdiabatic11.*Awave2;
    
    %Dwave1 = Dwave1.*exp(-XX.^2/50000-YY.^2/50000);
    %Dwave2 = Dwave2.*exp(-XX.^2/50000-YY.^2/50000);
    
    %Dwave1 = Dwave1.*(.5-(XX/LL/2).^2-(YY/LL/2).^2)*2;
    %exp(-(XX/50).^2-(YY/50).^2);
    %Dwave2 = Dwave2.*(.5-(XX/LL/2).^2-(YY/LL/2).^2)*2;
    %exp(-(XX/50).^2-(YY/50).^2);
    
    conservation(ll) = sum(sum(Dwave2.*conj(Dwave2)))+sum(sum(Dwave1.*conj(Dwave1)));
    
    %Dwave1 = Dwave1.*exp(-(XX/100).^2-(YY/100).^2);
    %Dwave2 = Dwave2.*exp(-(XX/100).^2-(YY/100).^2);
    
    %Dwave1 = Dwave1/conservation(ll);
    %Dwave2 = Dwave2/conservation(ll);
    
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
save('C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\G-k1-l_001-n_01-eV4-m1-LL50-NN8-dt5000.mat', '-v7.3')
%save(['C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])
%save(['C:\Chelsea\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])

x=linspace(0,2001,200);
y=sin(x);
sound(y)
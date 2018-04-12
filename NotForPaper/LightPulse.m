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

LL = 10;
NN = 2^8;
xx = linspace(-LL, LL-2*LL/NN, NN); %Angstrom
dt = 1/100000; %femtoseconds
totalTime = 0.1; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[XX,YY] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

%% Constants

kappa = 20; %eV/Angstrom
lambda = 10; %eV/Angstrom
mm = 10*0.0096; %amu
hbar = 0.658;
lighteV = 10; %eV
nu = 0; %eV/Angstrom
xposition = 5;
width = .5;
carrierphase = 0;
pulsetime = 0.3;
pulsewidth = 0.2;

omegaLight = lighteV/hbar;

%% Hamiltonians

W11 = kappa*XX;
%W11 = kappa/2*(XX.^2);
%.*(1-tanh(XX-50))
%W11 = kappa*(XX)+YY.^2+XX.^2;
%W11 = 0.5*(-kappa*(5-XX).*(1-tanh(-5+XX))+(YY/5).^2+(XX/10).^2);
%W11 = kappa*(-.1*(XX-5).^2+2.5+YY.^2/10+XX.^2/10);
W12 = lambda*YY;
%W22 = kappa/2*(XX.^2+YY.^2)+1;
%W22 = -kappa*(XX)+YY.^2+XX.^2;
%W22 = 0.5*(-kappa*(XX-5).*(1-tanh(XX-5))+(YY/5).^2+(XX/10).^2)-1;
%W22 = kappa*(-XX+YY.^2/10+XX.^2/10);
W22 = -kappa*XX;

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

%{
A = zeros(NN,NN,2,2);
Ainv = A;
WD = A;

for nn = 1:NN
    for mm = 1:NN
        
        [A(nn,mm,:,:), WD(nn,mm,:,:), Ainv(nn,mm,:,:)] = eig([W11(nn,mm), W12(nn,mm); W12(nn,mm), W22(nn, mm)]);
                
    end
end

DiabaticToAdiabatic11 = ((YY==0&XX<0)*-2+1).*A(:,:,1,1);
DiabaticToAdiabatic12 = ((YY==0&XX<0)*-2+1).*abs(A(:,:,1,2));
DiabaticToAdiabatic22 = DiabaticToAdiabatic11;
DiabaticToAdiabatic21 = -DiabaticToAdiabatic12;

%DiabaticToAdiabatic11(YY~=0&XX<0) = -DiabaticToAdiabatic11(YY~=0&XX<0);

AdiabaticToDiabatic11 = DiabaticToAdiabatic11;
AdiabaticToDiabatic12 = DiabaticToAdiabatic21;
AdiabaticToDiabatic22 = DiabaticToAdiabatic22;
AdiabaticToDiabatic21 = DiabaticToAdiabatic12;


%}

DiabaticToAdiabatic22 = abs(cos(0.5*atan(W12./(0.5*(W11-W22)))));
DiabaticToAdiabatic21 = abs(sin(0.5*atan(W12./(0.5*(W11-W22)))));
DiabaticToAdiabatic22(XX==0&YY>=0) = 1/sqrt(2);
DiabaticToAdiabatic22(XX==0&YY<0) = -1/sqrt(2);
%DiabaticToAdiabatic22 = sign(YY).*DiabaticToAdiabatic22;
DiabaticToAdiabatic11 = DiabaticToAdiabatic22;
DiabaticToAdiabatic12= -DiabaticToAdiabatic21;

AdiabaticToDiabatic11 = DiabaticToAdiabatic11;
AdiabaticToDiabatic22 = DiabaticToAdiabatic22;
AdiabaticToDiabatic12 = DiabaticToAdiabatic21;
AdiabaticToDiabatic21 = DiabaticToAdiabatic12;
%}
%{
W = 0.5*(W11-W22);

DiabaticToAdiabatic22 = sqrt(1./sqrt(W12.^2./(W).^2+1)+1)/sqrt(2);
DiabaticToAdiabatic22 = sign(YY).*DiabaticToAdiabatic22;
DiabaticToAdiabatic22(XX==0&YY>=0) = 1/sqrt(2);
DiabaticToAdiabatic22(XX==0&YY<0) = -1/sqrt(2);
DiabaticToAdiabatic11 = DiabaticToAdiabatic22;

DiabaticToAdiabatic21 = abs(W12)./(sqrt(2)*abs(W).*sqrt(1./sqrt(W12.^2./W.^2+1)+1).*sqrt((W12.^2+W.^2)./W.^2));
DiabaticToAdiabatic21(XX==0) = 1/sqrt(2);
DiabaticToAdiabatic12 = -DiabaticToAdiabatic21;

AdiabaticToDiabatic22 = DiabaticToAdiabatic22;
AdiabaticToDiabatic11 = DiabaticToAdiabatic11;
AdiabaticToDiabatic21 = DiabaticToAdiabatic12;
AdiabaticToDiabatic12 = DiabaticToAdiabatic21;
%}


%% Initial wavefunctions

%Dwave1 = exp(-mm/2/hbar*((XX+2).^2+(YY).^2)/.1);
Awave2 = exp(-(((XX+xposition)/width).^2+(YY/width).^2));
Awave2 = Awave2/sqrt(sum(sum(Awave2.*Awave2)));
Awave1 = zeros(size(Awave2));
%Dwave2 = exp(-(((XX+xposition)/width).^2+(YY/width).^2));
%Dwave2 = Dwave2/sqrt(sum(sum(Dwave2.*Dwave2)));
%Dwave1 = zeros(size(Dwave2));

%D0 = Dwave2;

Dwave1 = AdiabaticToDiabatic11.*Awave1 + AdiabaticToDiabatic12.*Awave2;
Dwave2 = AdiabaticToDiabatic21.*Awave1 + AdiabaticToDiabatic22.*Awave2;

zlimit = max(max(abs(Dwave2)))^2.*1.1;
zlimit2 = max(sum(Dwave2.*conj(Dwave2),1))*1.1;

%% Stuff for figures

wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
wavemovie2(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
anglemovie = wavemovie;
wavemovieD = wavemovie; 
jj = 100;

%% Propagate

time = 0;
conservation = zeros(tNN,1);
Excited = conservation;
Ground = conservation;
Dipolet = conservation;

disp('Coming up with time estimate...')

tic

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,jj) == 1
        
        %{
        f = figure();
        subplot(2,1,1)
        surf(XX, YY, Awave2.*conj(Awave2))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['Adiabatic excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX, YY, Awave1.*conj(Awave1))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['Adiabatic ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovie((ll-1)/jj+1) = getframe(f);
        %}
        
        e = figure();
        subplot(2,1,1)
        surf(XX, YY, Dwave2.*conj(Dwave2))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['Diabatic excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX, YY, Dwave1.*conj(Dwave1))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['Diabatic ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovieD((ll-1)/jj+1) = getframe(e);
        %}
        
        %{
        g = figure();
        subplot(2,2,1)
        plot(XX(1,:), sum(Awave2.*conj(Awave2),1))
        axis([-LL LL 0 zlimit2])
        title(['Adiabatic excited X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,2)
        plot(YY(:,1), sum(Awave2.*conj(Awave2),2))
        axis([-LL LL 0 zlimit2])
        title(['Adiabatic excited Y sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,3)
        plot(XX(1,:), sum(Awave1.*conj(Awave1),1))
        axis([-LL LL 0 zlimit2])
        title(['Adiabatic ground X sums ' num2str(ll*dt) ' fs'])
        subplot(2,2,4)
        plot(YY(:,1), sum(Awave1.*conj(Awave1),2))
        axis([-LL LL 0 zlimit2])
        title(['Adiabatic ground Y sums ' num2str(ll*dt) ' fs'])
        wavemovie2((ll-1)/jj+1) = getframe(g);
                
        bool = abs(Awave2) > 0.0001;
        Angle2 = angle(Awave2);
        Angle2 = bool.*Angle2;
        
        bool = abs(Awave1) > 0.0001;
        Angle1 = angle(Awave1);
        Angle1 = bool.*Angle1;
        
        h = figure();
        subplot(2,1,1)
        pcolor(XX, YY, Angle2.*conj(Angle2))
        axis([-LL LL -LL LL])
        caxis([0,2*pi])
        xlabel('x')
        ylabel('y')
        title(['Adiabatic excited ' num2str(ll*dt) ' fs'])
        shading flat
        colorbar
        subplot(2,1,2)
        pcolor(XX, YY, Angle1.*conj(Angle1))
        axis([-LL LL -LL LL])
        caxis([0,2*pi])
        xlabel('x')
        ylabel('y')
        title(['Adiabatic ground ' num2str(ll*dt) ' fs'])
        shading flat
        colorbar
        anglemovie((ll-1)/jj+1) = getframe(h);
        %}
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
    
    Dipolet(ll) = nu*sin(omegaLight*((ll-1)*dt-carrierphase));
    %nu*tanh((ll-1)*dt-pulsetime)*exp(-(((ll-1)*dt-pulsetime)/pulsewidth)^2);
    %nu*exp(-(((ll-1)*dt-pulsetime)/pulsewidth)^2);
    %nu*sin(omegaLight*((ll-1)*dt-carrierphase));
    
    Dipole =  Dipolet(ll)+zeros(size(Dwave1));
        
    UDipole11 = cos(2*abs(Dipole)*dt/2/hbar);
    UDipole12 = -1i*sin(2*abs(Dipole)*dt/2/hbar).*sign(Dipole);

    %% Propagation
    
    Dwave1 = ifft2(KEUU.*fft2(Dwave1));
    Dwave2 = ifft2(KEUU.*fft2(Dwave2));
    
    Dwave1 = U11.*Dwave1 + U12.*Dwave2;
    Dwave2 = U22.*Dwave2 + U12.*Dwave1;
    %{
    Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic12.*Dwave2;
    Awave2 = DiabaticToAdiabatic21.*Dwave1 + DiabaticToAdiabatic22.*Dwave2;
    
    Awave1 = UDipole11.*Awave1 + UDipole12.*Awave2;
    Awave2 = UDipole11.*Awave2 + UDipole12.*Awave1;
    
    Dwave1 = AdiabaticToDiabatic11.*Awave1 + AdiabaticToDiabatic12.*Awave2;
    Dwave2 = AdiabaticToDiabatic21.*Awave1 + AdiabaticToDiabatic22.*Awave2;
    %}
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
    
    %Excited(ll) = sum(sum(Awave2.*conj(Awave2)));
    %Ground(ll) = sum(sum(Awave1.*conj(Awave1)));
    Excited(ll) = sum(sum(Dwave2.*conj(Dwave2)));
    Ground(ll) = sum(sum(Dwave1.*conj(Dwave1)));
    
    
    
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

%save('C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\G-k1-l_001-n_01-eV4-m1-LL50-NN8-dt5000.mat', '-v7.3')

x=linspace(0,2001,200);
y=sin(x);
sound(y)
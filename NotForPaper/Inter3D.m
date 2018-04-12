%This one is happy to take light

clear all
close all
clc
disp('Coming up with time estimate...')

%% Discritize the wavefunction

LL = 10; %Angstrom
NN = 2^8; 
xx = linspace(-LL, LL, NN); %Angstrom
dt = 1/10; %femtoseconds
totalTime = 25; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

kk = linspace(-1/dx, 1/dx, NN);
kk = fftshift(kk);

[XX,YY,ZZ] = meshgrid(xx,xx,xx);
[k1,k2,k3] = meshgrid(kk,kk,kk);

%% Constants

kappax1 = .2;
%.25; %eV/Angstrom
kappax2 = .25;
%.0556;
kappay1 = .1; %eV/Angstrom
kappay2 = .1;
kappaz1 = .1;
kappaz2 = .1;
lambda = 0.01;
mm = 1*104.4070; %amu
%mm = 1/100;
hbar = 0.658;
lighteV = 0.5; %eV
nu = 0.1; %eV
xposition = -1; %Angstrom
width = 1;
%1/sqrt(kappa*mm)/2/hbar %Angstrom
carrierphase = 0; 
pulsedelay = 0; %fs
Woffset = 5;
XX01 = -1.0557;
%0;
XX02 = .944;
%2

omegax1 = (kappax1/mm)^.5;
widthx1 = 2*hbar/mm/omegax1;
omegax2 = (kappax2/mm)^.5;
widthx2 = 2*hbar/mm/omegax2;
omegay1 = (kappay1/mm)^.5;
widthy1 = 2*hbar/mm/omegay1;
omegay2 = (kappay2/mm)^.5;
widthy2 = 2*hbar/mm/omegay2;
omegaz1 = (kappaz1/mm)^.5;
widthz1 = 2*hbar/mm/omegaz1;
omegaz2 = (kappaz2/mm)^.5;
widthz2 = 2*hbar/mm/omegaz2;

E=.5;

omegaLight = lighteV/hbar;

%% Hamiltonians

W22 = kappax2*(XX-XX02).^2+kappay2*YY.^2+kappaz2*ZZ.^2-.223;
%W22 = .5*kappax2*XX.^2+.5*kappay2*YY.^2+.5*kappaz2*ZZ.^2;
%W11 = .5*kappax1*XX.^2+.5*kappay1*YY.^2+.5*kappaz1*ZZ.^2;
%W11 = kappax1*(-XX+XX01).^2+kappay1*YY.^2+kappaz1*ZZ.^2-1;
W11 = kappax1*(XX-XX01).^2+kappay1*YY.^2+kappaz1*ZZ.^2-.223;
W12 = lambda*YY;
%W22 = .5+zeros(size(XX));
%W11 = -.5+zeros(size(XX));
%W12 = zeros(size(XX));


KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2+k3.^2));

%KEUU1 = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2));
%KEUU2 = exp(-1i*pi^2*hbar/8/mm*dt*(k2.^2));
%KEUU3 = exp(-1i*pi^2*hbar/8/mm*dt*(k3.^2));

C = (W11-W22).^2;
A = exp(-1i/hbar*(W11+W22)*dt/2);
B = (W22-W11);

D = sqrt(4*W12.*conj(W12)+C);

U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
U12 = A.*1i.*sin(D*dt/2/hbar)./D*(-2).*W12;

%% Initial wavefunctions

Dwave2 = exp(-((XX-xposition).^2/widthx2+YY.^2/widthy2+ZZ.^2/widthz2));

%Dwave2 = exp(-sqrt(kappa*mm)/2/hbar*(XX.^2+YY.^2));

Dwave2 = Dwave2/sqrt(sum(sum(sum(Dwave2.*Dwave2))));
Dwave1 = zeros(size(Dwave2));

D0 = Dwave2;

zlimit = max(max(max(abs(Dwave2))))^2.*1.5;
zlimit2 = max(sum(sum(Dwave2.*conj(Dwave2),1)))*1.5;

%% Stuff for figures

wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
wavemovie2 = wavemovie;
anglemovie = wavemovie;
wavemovieA = wavemovie;
wavemovie2A = wavemovie;
anglemovieA = wavemovie;
jj = 10;

%% Initiate

time = 0;
conservation = zeros(tNN,1);
Excited = conservation;
Ground = conservation;
Dipolet = conservation;
Xcenter = conservation;
Ycenter = conservation;
Zcenter = conservation;
PE = conservation;
KE = conservation;
TE = conservation;

%% Adiabatic to diabatic
AA = zeros(NN,NN,NN,2,2);
Ainv = AA;
WD = AA;

for nn = 1:NN
    for pp = 1:NN
        for ll = 1:NN
        
            [AA(nn,pp,ll,:,:), WD(nn,pp,ll,:,:), Ainv(nn,pp,ll,:,:)] = eig([W11(nn,pp,ll), W12(nn,pp,ll); W12(nn,pp,ll), W22(nn,pp,ll)]);
        end     
    end
end

DiabaticToAdiabatic11 = AA(:,:,:,1,1);
DiabaticToAdiabatic12 = AA(:,:,:,1,2);
DiabaticToAdiabatic22 = AA(:,:,:,2,2);
DiabaticToAdiabatic21 = AA(:,:,:,2,1);

AdiabaticToDiabatic11 = Ainv(:,:,:,1,1);
AdiabaticToDiabatic12 = Ainv(:,:,:,1,2);
AdiabaticToDiabatic21 = Ainv(:,:,:,2,1);
AdiabaticToDiabatic22 = Ainv(:,:,:,2,2);

%Dwave1 = AdiabaticToDiabatic11.*Awave1 + AdiabaticToDiabatic12.*Awave2;
%Dwave2 = AdiabaticToDiabatic21.*Awave1 + AdiabaticToDiabatic22.*Awave2;

Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic12.*Dwave2;
Awave2 = DiabaticToAdiabatic21.*Dwave1 + DiabaticToAdiabatic22.*Dwave2;

tic
%% Propagation

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,jj) == 1
        
        
        f = figure();
        subplot(2,1,1)
        surf(XX(:,:,NN/2), YY(:,:,NN/2), Dwave2(:,:,NN/2).*conj(Dwave2(:,:,NN/2)))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX(:,:,NN/2), YY(:,:,NN/2), Dwave1(:,:,NN/2).*conj(Dwave1(:,:,NN/2)))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovie((ll-1)/jj+1) = getframe(f);
        
        
        
        F = figure();
        subplot(2,1,1)
        surf(XX(:,:,NN/2), YY(:,:,NN/2), Awave2(:,:,NN/2).*conj(Awave2(:,:,NN/2)))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['A excited ' num2str(ll*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX(:,:,NN/2), YY(:,:,NN/2), Awave1(:,:,NN/2).*conj(Awave1(:,:,NN/2)))
        axis([-LL LL -LL LL 0 zlimit])
        xlabel('x')
        ylabel('y')
        title(['A ground ' num2str(ll*dt) ' fs'])
        shading flat
        wavemovieA((ll-1)/jj+1) = getframe(F);
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
    
    if (ll-1)*dt >= pulsedelay
        Dipolet(ll) = nu*sin(omegaLight*(((ll-1)*dt-pulsedelay)-carrierphase));
    else
        Dipolet(ll) = 0;
    end
    
    UDipole11 = cos(2*abs(Dipolet(ll))*dt/2/hbar);
    UDipole12 = -1i*sin(2*abs(Dipolet(ll))*dt/2/hbar).*sign(Dipolet(ll));
        
    %D = sqrt(4*W12.*conj(W12)+C);

    U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
    U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
    U12 = A.*1i.*sin(D*dt/2/hbar)./D.*(-2).*W12;

    %% Propagation
    
    FFTD2 = fftn(Dwave2);
    FFTD1 = fftn(Dwave1);
    
    Dwave1 = ifftn(KEUU.*FFTD1);
    Dwave2 = ifftn(KEUU.*FFTD2);
    %{
    Dwave1 = ifft(KEUU1.*fft(Dwave1,[],1),[],1);
    Dwave1 = ifft(KEUU2.*fft(Dwave1,[],2),[],2);
    Dwave1 = ifft(KEUU3.*fft(Dwave1,[],3),[],3);
    
    Dwave2 = ifft(KEUU1.*fft(Dwave2,[],1),[],1);
    Dwave2 = ifft(KEUU2.*fft(Dwave2,[],2),[],2);
    Dwave2 = ifft(KEUU3.*fft(Dwave2,[],3),[],3);
    %}
    Dwave1 = U11.*Dwave1 + U12.*Dwave2;
    Dwave2 = U22.*Dwave2 + U12.*Dwave1;
    
    Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic12.*Dwave2;
    Awave2 = DiabaticToAdiabatic21.*Dwave1 + DiabaticToAdiabatic22.*Dwave2;
    
    Awave1 = UDipole11.*Awave1 + UDipole12.*Awave2;
    Awave2 = UDipole11.*Awave2 + UDipole12.*Awave1;
    
    Dwave1 = AdiabaticToDiabatic11.*Awave1 + AdiabaticToDiabatic12.*Awave2;
    Dwave2 = AdiabaticToDiabatic21.*Awave1 + AdiabaticToDiabatic22.*Awave2;
    %}
    conservation(ll) = sum(sum(sum(Dwave2.*conj(Dwave2)))+sum(sum(Dwave1.*conj(Dwave1))));
    Excited(ll) = sum(sum(sum(Awave2.*conj(Awave2))));
    Ground(ll) = sum(sum(sum(Awave1.*conj(Awave1))));
    PE(ll) = real(sum(sum(sum(Dwave2.*W22.*conj(Dwave2)))));
    KE(ll) = real(sum(sum(sum(pi^2*hbar^2/8/mm*FFTD2.*(k1.^2+k2.^2+k3.^2).*conj(FFTD2)))/NN^3));
    TE(ll) = PE(ll)+KE(ll);
        
    Xcenter(ll) = real(sum(sum(sum(XX.*Dwave2.*conj(Dwave2)))));
    Ycenter(ll) = real(sum(sum(sum(YY.*Dwave2.*conj(Dwave2)))));
    Zcenter(ll) = real(sum(sum(sum(ZZ.*Dwave2.*conj(Dwave2)))));
    
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

%save('C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\TestInter10fs.mat', '-v7.3')

x=linspace(0,2001,200);
y=sin(x);
sound(y)
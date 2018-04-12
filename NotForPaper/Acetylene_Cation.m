EV = linspace(2,1,3);
DD = length(EV);

load('CationPES.mat')

%for TT = 1:DD
TT=1;
%clearvars -except EV TT DD
close all
clc
disp('Coming up with time estimate...')
disp(['EV ' num2str(EV(TT))])

%% Discritize the wavefunction

LL = 3; %Angstrom
NN = 2^9; 
xx = linspace(-LL, LL, NN); %Angstrom
stretch = linspace(0.8,20,NN);
dt = 1/10; %femtoseconds
totalTime = 15; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);
%tNN=1;

kk = linspace(-1/dx, 1/dx, NN);
kk = fftshift(kk);

%[XX,YY] = meshgrid(stretch,xx);
[XX,YY] = meshgrid(stretch,xx);
[k1,k2] = meshgrid(kk,kk);

%% Constants

kappax1 = 8;
kappax2 = 8;
kappay1 = 8; %eV/Angstrom
kappay2 = 8;
lambda = .2;
mm = 1*104.4070; %amu
hbar = 0.658;
lighteV = EV(TT); %eV
mu = 0.5*(XX<3).*(YY<3); %eV
xposition = 2; %Angstrom
carrierphase = 0; 
pulsedelay = 0; %fs

omegax1 = (kappax1/mm)^.5;
widthx1 = 2*hbar/mm/omegax1;
omegax2 = (kappax2/mm)^.5;
widthx2 = 2*hbar/mm/omegax2;
omegay1 = (kappay1/mm)^.5;
widthy1 = 2*hbar/mm/omegay1;
omegay2 = (kappay2/mm)^.5;
widthy2 = 2*hbar/mm/omegay2;

omegaLight = lighteV/hbar;

%% Hamiltonians
CC_ASig = CC_ASig*0.529177249;
CC_XPi = CC_XPi*0.529177249;
[c1,a1,b1] = unique(CC_ASig);
SigNew = smooth(CC_ASig,ASigma,30);
W22 = spline(CC_ASig(a1), SigNew(a1), XX)+kappay2*YY.^2;
[c1,a1,b1] = unique(CC_XPi);
XPiNew = smooth(CC_XPi,XPi,30);
W11 = spline(CC_XPi(a1), XPiNew(a1), XX)+kappay1*YY.^2;
W12 = lambda*YY;

KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2));

C = (W11-W22).^2;
A = exp(-1i/hbar*(W11+W22)*dt/2);
B = (W22-W11);

D = sqrt(4*W12.*conj(W12)+C);

U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
U12 = A.*1i.*sin(D*dt/2/hbar)./D*(-2).*W12;

%% Initial wavefunctions

Dwave2 = exp(-((XX-xposition).^2/widthx2+YY.^2/widthy2));
Dwave2 = Dwave2/sqrt(sum(sum(sum(Dwave2.*Dwave2))));
Dwave1 = zeros(size(Dwave2));

D0 = Dwave2;

zlimit = max(max(max(abs(Dwave2))))^2.*1.5;
zlimit2 = max(sum(sum(Dwave2.*conj(Dwave2),1)))*1.5;

%% Stuff for figures

jj = 5;

wavemovie(1:(floor(tNN/jj)+1)) = struct('cdata',[],'colormap',[]);
wavemovie2 = wavemovie;
anglemovie = wavemovie;
wavemovieA = wavemovie;
wavemovie2A = wavemovie;
anglemovieA = wavemovie;


%% Initiate

time = 0;
conservation = zeros(tNN,1);
ExcitedA = conservation;
GroundA = conservation;
ExcitedD = conservation;
GroundD = conservation;
XcenterA1 = conservation;
XcenterA2 = conservation;
XcenterD1 = conservation;
XcenterD2 = conservation;
PE = conservation;
KE = conservation;
TE = conservation;

%% Adiabatic to diabatic
AA = zeros(NN,NN,2,2);
Ainv = AA;
WD = AA;

for nn = 1:NN
    for pp = 1:NN

        [AA(nn,pp,:,:), WD(nn,pp,:,:)] = eig([W11(nn,pp), W12(nn,pp); W12(nn,pp), W22(nn,pp)]);
        Ainv(nn,pp,:,:) = inv(squeeze(AA(nn,pp,:,:)));

    end
end

DiabaticToAdiabatic11 = AA(:,:,1,1);
DiabaticToAdiabatic12 = AA(:,:,1,2);
DiabaticToAdiabatic22 = AA(:,:,2,2);
DiabaticToAdiabatic21 = AA(:,:,2,1);

AdiabaticToDiabatic11 = Ainv(:,:,1,1);
AdiabaticToDiabatic12 = Ainv(:,:,1,2);
AdiabaticToDiabatic21 = Ainv(:,:,2,1);
AdiabaticToDiabatic22 = Ainv(:,:,2,2);

Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic12.*Dwave2;
Awave2 = DiabaticToAdiabatic21.*Dwave1 + DiabaticToAdiabatic22.*Dwave2;

tic
%% Propagation

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,jj) == 1
        
        f = figure();
        subplot(2,1,1)
        pcolor(XX, YY, Dwave2.*conj(Dwave2))
        xlabel('x')
        ylabel('y')
        title(['excited ' num2str(ll*dt) ' fs'])
        shading flat
        caxis([0 zlimit])
        axis([0.8,20,-3,3])
        subplot(2,1,2)
        pcolor(XX, YY, Dwave1.*conj(Dwave1))
        xlabel('x')
        ylabel('y')
        title(['ground ' num2str(ll*dt) ' fs'])
        shading flat
        caxis([0 zlimit])
        axis([0.8,20,-3,3])
        set(f, 'position', [100,50,450,850])
        wavemovie((ll-1)/jj+1) = getframe(gcf);
        
        F = figure();
        subplot(2,1,1)
        pcolor(XX, YY, Awave2.*conj(Awave2))
        xlabel('x')
        ylabel('y')
        title(['A excited ' num2str(ll*dt) ' fs'])
        shading flat
        caxis([0 zlimit*.05])
        axis([0.8,20,-3,3])
        subplot(2,1,2)
        pcolor(XX, YY, Awave1.*conj(Awave1))
        xlabel('x')
        ylabel('y')
        title(['A ground ' num2str(ll*dt) ' fs'])
        shading flat
        caxis([0 zlimit*.05])
        axis([0.8,20,-3,3])
        set(F, 'position', [100,50,450,850])
        wavemovieA((ll-1)/jj+1) = getframe(gcf);
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
    
    if omegaLight ~=0
        Dipolet = mu*sin(omegaLight*(((ll-1)*dt-pulsedelay)-carrierphase));
    elseif omegaLight == 0
        Dipolet = mu;
    end
    
    UDipole11 = cos(2*abs(Dipolet)*dt/2/hbar);
    UDipole12 = -1i*sin(2*abs(Dipolet)*dt/2/hbar).*sign(Dipolet);

    %% Propagation
    
    FFTD2 = fftn(Dwave2);
    FFTD1 = fftn(Dwave1);
    
    Dwave1 = ifftn(KEUU.*FFTD1);
    Dwave2 = ifftn(KEUU.*FFTD2);

    Dwave1 = U11.*Dwave1 + U12.*Dwave2;
    Dwave2 = U22.*Dwave2 + U12.*Dwave1;
    
    Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic12.*Dwave2;
    Awave2 = DiabaticToAdiabatic21.*Dwave1 + DiabaticToAdiabatic22.*Dwave2;
    
    Awave1 = UDipole11.*Awave1 + UDipole12.*Awave2;
    Awave2 = UDipole11.*Awave2 + UDipole12.*Awave1;
    
    Dwave1 = AdiabaticToDiabatic11.*Awave1 + AdiabaticToDiabatic12.*Awave2;
    Dwave2 = AdiabaticToDiabatic21.*Awave1 + AdiabaticToDiabatic22.*Awave2;

    conservation(ll) = sum(sum(sum(Dwave2.*conj(Dwave2)))+sum(sum(Dwave1.*conj(Dwave1))));
    ExcitedA(ll) = sum(sum(sum(Awave2.*conj(Awave2))));
    GroundA(ll) = sum(sum(sum(Awave1.*conj(Awave1))));
    ExcitedD(ll) = sum(sum(sum(Dwave2.*conj(Dwave2))));
    GroundD(ll) = sum(sum(sum(Dwave1.*conj(Dwave1))));
    PE(ll) = real(sum(sum(sum(Dwave2.*W22.*conj(Dwave2)))));
    KE(ll) = real(sum(sum(sum(pi^2*hbar^2/8/mm*FFTD2.*(k1.^2+k2.^2).*conj(FFTD2)))/NN^3));
    TE(ll) = PE(ll)+KE(ll);
        
    XcenterA2(ll) = real(sum(sum(sum(Awave2.*XX.*conj(Awave2)))));
    XcenterA1(ll) = real(sum(sum(sum(Awave1.*XX.*conj(Awave1)))));
    XcenterD1(ll) = real(sum(sum(sum(Dwave1.*XX.*conj(Dwave1)))));
    XcenterD2(ll) = real(sum(sum(sum(Dwave2.*XX.*conj(Dwave2)))));
    
    if ll == 1
    
        seconds = (tNN-ll)*(toc+(1-jj)/jj*second);
        hours = floor(seconds/60/60);
        minutes = floor((seconds - hours*60*60)/60);
        seconds = seconds - hours*60*60 - minutes*60;
        disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])
        
    end
end

%% Done :)

t = datetime;

%save(['C:\Chelsea\TDSE_Output\BigFix\NewFigures_EV' num2str(EV(TT)) '_Pos-1.25_Mu0.1.mat'], '-v7.3')

%end

x=linspace(0,2001,200);
y=sin(x);
sound(y)
%}
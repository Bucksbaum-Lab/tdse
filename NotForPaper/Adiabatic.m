clearvars

Position = -1.25;
Nu = 0.05;
EV = 0.5;

close all
clc
disp('Coming up with time estimate...')

%% Discritize the wavefunction

LL = 10; %Angstrom
NN = 2^8; 
xx = linspace(-LL, LL, NN); %Angstrom
dt = 1/10; %femtoseconds
totalTime = 30; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

kk = linspace(-1/dx, 1/dx, NN);
kk = fftshift(kk);

[XX,YY,ZZ] = meshgrid(xx,xx,xx);
[k1,k2,k3] = meshgrid(kk,kk,kk);

%% Constants

kappax1 = .2; %eV/Angstrom
kappax2 = .25;
kappay1 = .1; %eV/Angstrom
kappay2 = .1;
kappaz1 = .1;
kappaz2 = .1;
lambda = 0.01;
mm = 1*104.4070; %amu
hbar = 0.658;
lighteV = EV; %eV
nu = Nu.*(XX<3).*(YY<3).*(ZZ<3); %eV
xposition = Position; %Angstrom
carrierphase = 0; 
pulsedelay = 0; %fs
Woffset = 5;
XX01 = -1.0557;
XX02 = .944;

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

omegaLight = lighteV/hbar;

%% Hamiltonians

W22 = kappax2*(XX-XX02).^2+kappay2*YY.^2+kappaz2*ZZ.^2-.223;
W11 = kappax1*(XX-XX01).^2+kappay1*YY.^2+kappaz1*ZZ.^2-.223;
W12 = lambda*YY;

KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2+k3.^2));

%% Initial wavefunctions

Dwave2 = exp(-((XX-xposition).^2/widthx2+YY.^2/widthy2+ZZ.^2/widthz2));
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
        
            [AA(nn,pp,ll,:,:), WD(nn,pp,ll,:,:)] = eig([W11(nn,pp,ll), W12(nn,pp,ll); W12(nn,pp,ll), W22(nn,pp,ll)]);
        end     
    end
end

DiabaticToAdiabatic11 = AA(:,:,:,1,1);
DiabaticToAdiabatic12 = AA(:,:,:,1,2);
DiabaticToAdiabatic22 = AA(:,:,:,2,2);
DiabaticToAdiabatic21 = AA(:,:,:,2,1);

Awave1 = DiabaticToAdiabatic11.*Dwave1 + DiabaticToAdiabatic12.*Dwave2;
Awave2 = DiabaticToAdiabatic21.*Dwave1 + DiabaticToAdiabatic22.*Dwave2;

Vp = WD(:,:,:,1,1);
Vm = WD(:,:,:,2,2);

A = exp(-1i/hbar*(Vp+Vm)*dt/2);
B = (Vm-Vp);
D = abs(Vp-Vm);

U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
U12 = A.*1i.*sin(D*dt/2/hbar)./D*(-2).*W12;

tic
%% Propagation

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,jj) == 1
        
        F = figure();
        subplot(2,1,1)
        pcolor(XX(:,:,NN/2), YY(:,:,NN/2), Awave2(:,:,NN/2).*conj(Awave2(:,:,NN/2)))
        xlabel('x')
        ylabel('y')
        title(['A excited ' num2str(ll*dt) ' fs'])
        shading flat
        caxis([0 zlimit*.05])
        axis([-3,3,-3,3])
        subplot(2,1,2)
        pcolor(XX(:,:,NN/2), YY(:,:,NN/2), Awave1(:,:,NN/2).*conj(Awave1(:,:,NN/2)))
        xlabel('x')
        ylabel('y')
        title(['A ground ' num2str(ll*dt) ' fs'])
        shading flat
        caxis([0 zlimit*.05])
        axis([-3,3,-3,3])
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
        Dipolet = nu*sin(omegaLight*(((ll-1)*dt-pulsedelay)-carrierphase));
    elseif omegaLight == 0
        Dipolet = nu;
    end
    
    UDipole11 = cos(2*abs(Dipolet)*dt/2/hbar);
    UDipole12 = -1i*sin(2*abs(Dipolet)*dt/2/hbar).*sign(Dipolet);
    
    U11 = A.*(cos(D*dt/2/hbar)+1i*sin(D*dt/2/hbar)./D.*B);
    U22 = A.*(cos(D*dt/2/hbar)-1i*sin(D*dt/2/hbar)./D.*B);
    U12 = zeros(size(U11));

    %% Propagation
    
    FFTA2 = fftn(Awave2);
    FFTA1 = fftn(Awave1);
    
    Awave1 = ifftn(KEUU.*FFTA1);
    Awave2 = ifftn(KEUU.*FFTA2);
   
    Awave1 = U11.*Awave1 + U12.*Awave2;
    Awave2 = U22.*Awave2 + U12.*Awave1;

    Awave1 = UDipole11.*Awave1 + UDipole12.*Awave2;
    Awave2 = UDipole11.*Awave2 + UDipole12.*Awave1;

    conservation(ll) = sum(sum(sum(Awave2.*conj(Awave2)))+sum(sum(Awave1.*conj(Awave1))));
    Excited(ll) = sum(sum(sum(Awave2.*conj(Awave2))));
    Ground(ll) = sum(sum(sum(Awave1.*conj(Awave1))));
    PE(ll) = real(sum(sum(sum(Awave2.*W22.*conj(Awave2)))));
    KE(ll) = real(sum(sum(sum(pi^2*hbar^2/8/mm*FFTA2.*(k1.^2+k2.^2+k3.^2).*conj(FFTA2)))/NN^3));
    TE(ll) = PE(ll)+KE(ll);
        
    Xcenter(ll) = real(sum(sum(sum(XX.*Awave2.*conj(Awave2)))));
    Ycenter(ll) = real(sum(sum(sum(YY.*Awave2.*conj(Awave2)))));
    Zcenter(ll) = real(sum(sum(sum(ZZ.*Awave2.*conj(Awave2)))));
    
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

save(['C:\Chelsea\TDSE_Output\BigFix\Adiabatic_EV' num2str(EV(TT)) '_Pos' num2str(Position(TT)) '_Nu' num2str(Nu(TT)) '.mat'], '-v7.3')

x=linspace(0,2001,200);
y=sin(x);
sound(y)
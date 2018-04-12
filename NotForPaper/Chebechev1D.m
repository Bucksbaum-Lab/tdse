clear all
close all
clc

%% Discritize the wavefunction

LL = 100;
NN = 2^10;
xx = linspace(-LL, LL-2*LL/NN, NN); %Angstrom
dt = 1/100; %femtoseconds
totalTime = 10; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

%% Constants

kappa = 1; %eV/Angstrom
lambda = 0; %eV/Angstrom
mm = 10*0.0096; %amu
hbar = 0.6582;
lighteV = .1; %eV
nu = 1; %eV/Angstrom

omegaLight = lighteV/hbar;

%% Hamiltonians

%W11 = kappa/2*xx.*xx;
%W11 = kappa*xx.^2/2;
W11 = zeros(size(xx));
W12 = lambda*xx;
%W22 = kappa*xx.^2+10;
W22 = W11+1;

%% Initial wavefunctions

Dwave1 = exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%+2*sqrt(sqrt(kappa*mm)/hbar)*xx.*exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%Dwave1 = sin(xx);
Dwave1 = Dwave1/sqrt(sum(sum(Dwave1.*Dwave1)));
%Dwave2 = Dwave1;
Dwave2 = zeros(size(Dwave1));
%Dwave2 = exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%Dwave2 = Dwave2/sqrt(sum(sum(Dwave2.*Dwave2)));
%Dwave1 = zeros(size(Dwave2));

D0 = Dwave1;

zlimit = (max(Dwave1.*conj(Dwave1))+ max(Dwave2.*conj(Dwave2)))*1.1;

ll = 25;
wavemovie(floor(tNN/ll)) = struct('cdata',[],'colormap',[]);

%% Propagate

nn = 100;
tic

conservation = zeros(1,tNN);
Difference1 = zeros(tNN, NN);
Difference2 = Difference1;

for ii = 1:tNN

    if mod(ii,ll) == 1
                
        f = figure();
        subplot(2,1,1)
        plot(xx, real(Dwave2.*conj(Dwave2)))
        axis([-LL LL-2*LL/NN 0 zlimit])
        title(['excited ' num2str(ii*dt) ' fs'])
        subplot(2,1,2)
        plot(xx, real(Dwave1.*conj(Dwave1)))
        axis([-LL LL-2*LL/NN 0 zlimit])
        title(['ground ' num2str(ii*dt) ' fs'])
        wavemovie((ii-1)/ll+1) = getframe(f);
        close all
        
        seconds = (tNN-ii)/ll*toc;
        tic
        hours = floor(seconds/60/60);
        minutes = floor((seconds - hours*60*60)/60);
        seconds = seconds - hours*60*60 - minutes*60;
        disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])
                
    end
    
    if (ii*dt) < 5
        [phi1, phi2, Emin, deltaEgrid] = GetPhi(W11, W22, W12+nu*sin(omegaLight*dt*(ii-1)), Dwave1, Dwave2, dx, mm, nn, xx);
    else
        [phi1, phi2, Emin, deltaEgrid] = GetPhi(W11, W22, W12, Dwave1, Dwave2, dx, mm, nn, xx);
    end
    A = GetA(deltaEgrid*dt/2/hbar, nn);
    
    Sum1 = zeros(1,size(phi1,2));
    Sum2 = Sum1;
    
    for jj = 1:nn
        Sum1 = A(jj)*phi1(jj,:)+Sum1;
        Sum2 = A(jj)*phi2(jj,:)+Sum2;
    end
        
    Dwave1 = exp(-1i*(Emin+deltaEgrid/2)*dt/hbar).*Sum1;
    Dwave2 = exp(-1i*(Emin+deltaEgrid/2)*dt/hbar).*Sum2;
    %-1i*1000*(sec(xx/10001*pi)-1)
    %-1i*(xx/100).^2
    
    conservation(ii) = sum(conj(Dwave1).*Dwave1+conj(Dwave2).*Dwave2);
    
    Difference1(ii,:) = Dwave1 - Dwave1/conservation(ii);
    Difference2(ii,:) = Dwave2 - Dwave2/conservation(ii);
    
    Dwave1 = Dwave1/conservation(ii);
    Dwave2 = Dwave2/conservation(ii);
    
end

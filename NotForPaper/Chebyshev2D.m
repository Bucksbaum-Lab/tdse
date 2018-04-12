clear all
close all
clc

%% Discritize the wavefunction

LL = 50;
NN = 2^8;
xx = linspace(-LL, LL-2*LL/NN, NN); %Angstrom
[XX,YY] = meshgrid(xx,xx);
dt = 1/100; %femtoseconds
totalTime = 10; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

%% Constants

kappa = .1; %eV/Angstrom
lambda = 0; %eV/Angstrom
mm = 10*0.0096; %amu
hbar = 0.6582;
lighteV = 1; %eV
nu = 0; %eV/Angstrom

omegaLight = lighteV/hbar;

%% Hamiltonians

%W11 = kappa/2*xx.*xx;
%W11 = kappa*XX;
W11 = kappa/2*(XX.^2+YY.^2);
W12 = lambda*YY;
W22 = -kappa*XX;
%W22 = zeros(size(W11));
%W22 = kappa/2*(XX.^2+YY.^2)+1;

%% Initial wavefunctions

Dwave1 = exp(-sqrt(kappa*mm)/2/hbar*((XX).^2+YY.^2));
%Dwave1 = exp(-sqrt(kappa*mm)/2/hbar*(xx.^2))+2*sqrt(sqrt(kappa*mm)/hbar)*xx.*exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%Dwave1 = sin(xx);
Dwave1 = Dwave1/sqrt(sum(sum(Dwave1.*Dwave1)));
Dwave2 = zeros(size(Dwave1));
%Dwave2 = exp(-sqrt(kappa*mm)/2/hbar*(xx.^2));
%Dwave2 = Dwave2/sqrt(sum(sum(Dwave2.*Dwave2)));
%Dwave1 = zeros(size(Dwave2));
D0 = Dwave1;

zlimit = max(max(abs(Dwave1)))^2.*1.1;
zlimit2 = max([sum(abs(Dwave1),1),(sum(abs(Dwave1),2))'])*1.1;

ll = 25;
wavemovie(floor(tNN/ll)) = struct('cdata',[],'colormap',[]);
wavemovie2 = wavemovie;

%% Propagate

nn = 33;
tic

conservation = zeros(1,tNN);
Energy = conservation;

disp('Coming up with time estimate...')

for ii = 1:tNN

    if mod(ii,ll) == 1
                
        f = figure();
        subplot(2,1,1)
        surf(XX, YY, Dwave2.*conj(Dwave2))
        axis([-50 50 -50 50 0 zlimit])
        title(['excited ' num2str(ii*dt) ' fs'])
        shading flat
        subplot(2,1,2)
        surf(XX, YY, Dwave1.*conj(Dwave1))
        axis([-50 50 -50 50 0 zlimit])
        title(['ground ' num2str(ii*dt) ' fs'])
        shading flat
        wavemovie((ii-1)/ll+1) = getframe(f);
        
        g = figure();
        subplot(2,2,1)
        plot(XX(NN/2,:), Dwave1(NN/2,:).*conj(Dwave1(NN/2,:)))
        axis([-50 50 0 zlimit])
        title(['ground X lineout ' num2str(ii*dt) ' fs'])
        subplot(2,2,2)
        plot(YY(:,NN/2), Dwave1(:,NN/2).*conj(Dwave1(:,NN/2)))
        axis([-50 50 0 zlimit])
        title(['ground Y lineout ' num2str(ii*dt) ' fs'])
        subplot(2,2,3)
        plot(XX(NN/2,:), Dwave2(NN/2,:).*conj(Dwave2(NN/2,:)))
        axis([-50 50 0 zlimit])
        title(['excited X lineout ' num2str(ii*dt) ' fs'])
        subplot(2,2,4)
        plot(YY(:,NN/2), Dwave2(:,NN/2).*conj(Dwave2(:,NN/2)))
        axis([-50 50 0 zlimit])
        title(['excited Y lineout ' num2str(ii*dt) ' fs'])
        wavemovie2((ii-1)/ll+1) = getframe(g);
        close all
        
        if ii ~= 1
            seconds = (tNN-ii)/ll*toc;
            tic
            hours = floor(seconds/60/60);
            minutes = floor((seconds - hours*60*60)/60);
            seconds = seconds - hours*60*60 - minutes*60;
            disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])
        else
            second = toc;
        end

    end
    
    [phi1, phi2, Emin, deltaEgrid] = GetPhi2D(W11, W22, W12+nu*sin(omegaLight*dt*(ii-1)), Dwave1, Dwave2, dx, mm, nn, xx);
    
    A = GetA(deltaEgrid*dt/2/hbar, nn);
    
    Sum1 = zeros(size(phi1,2), size(phi1,2),1);
    Sum2 = Sum1;
    
    for jj = 1:nn
        Sum1 = A(jj)*phi1(:,:,jj)+Sum1;
        Sum2 = A(jj)*phi2(:,:,jj)+Sum2;
    end
        
    Dwave1 = exp(-1i*(Emin+deltaEgrid/2)*dt/hbar).*Sum1;
    Dwave2 = exp(-1i*(Emin+deltaEgrid/2)*dt/hbar)*Sum2;
    
    conservation(ii) = sum(sum(conj(Dwave1).*Dwave1+conj(Dwave2).*Dwave2));
    Energy(ii) = sum(sum((-4*hbar^2/2/mm*del2(Dwave1,xx,xx)+W11.*Dwave1).*conj(Dwave1)));
    %{
    Dwave1(1,1) = 0;
    Dwave1(1,NN) = 0;
    Dwave1(NN,1) = 0;
    Dwave1(NN,NN) = 0;
    
    Dwave2(1,1) = 0;
    Dwave2(1,NN) = 0;
    Dwave2(NN,1) = 0;
    Dwave2 (NN,NN) = 0;
    %}
    %Dwave1 = Dwave1/conservation(ii);
    %Dwave2 = Dwave2/conservation(ii);
    
    if ii == 1
    
        seconds = (tNN-ii)*(toc+(1-ll)/ll*second);
        hours = floor(seconds/60/60);
        minutes = floor((seconds - hours*60*60)/60);
        seconds = seconds - hours*60*60 - minutes*60;
        disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])
        
    end
    
end

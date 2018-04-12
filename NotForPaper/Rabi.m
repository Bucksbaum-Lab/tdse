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
NN = 2^12;
xx = linspace(-LL, LL-2*LL/NN, NN); %Angstrom
dt = 1/100; %femtoseconds
totalTime = 10; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

%[XX,YY] = meshgrid(xx,xx);
%[k1,k2] = meshgrid(kk,kk);

%% Constants

kappa = 1; %eV/Angstrom^2
mm = 500*0.0096; %amu
hbar = 0.6582;  %eV*fs
lighteV = 1.2; %eV
nu = 1; %eV/Angstrom

omegaLight = lighteV/hbar;

%% Hamiltonians

tic
%W11 = 0.5*kappa*xx.^2;
W11 = (10*(xx).^2-(xx).^3)*0.01;
W12 = -nu*xx.*sin(omegaLight*0);
W22 = (20*(xx-1).^2)*0.01+1;
%W22 = (10*xx.^2-xx.^3-xx)/1000+10*.0096;

[V1, V2, A, Ainv] = makeW1D(W11, W12, W22);

%V1 = V1 + exp(-(XX.^2+YY.^2))*1i;
%V2 = V2 + exp(-(XX.^2+YY.^2))*1i;

KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(kk.^2+kk.^2)/2);

U1 = exp(-1i*V1/hbar*dt/2*1.512);
U2 = exp(-1i*V2/hbar*dt/2*1.512);

%% Initial wavefunctions

%Dwave1 = exp(-sqrt(kappa*mm)/(2*sqrt(2))/hbar*((xx).^2));
Dwave1 = exp(-((xx).^2));
Dwave1 = Dwave1/sqrt(sum(sum(Dwave1.*Dwave1)));
Dwave2 = zeros(size(Dwave1));
%Dwave1 = zeros(size(Dwave2));

Awave1 = (A(:,1,1))'.*Dwave1 + (A(:,1,2))'.*Dwave2;
Awave2 = (A(:,2,1))'.*Dwave1 + (A(:,2,2))'.*Dwave2;

zlimit = max(max(abs(Dwave1))).*1.1;
zlimit2 = max([sum(abs(Dwave1),1),(sum(abs(Dwave1),2))'])*1.1;

seconds = tNN*toc;
hours = floor(seconds/60/60);
minutes = floor((seconds - hours*60*60)/60);
seconds = seconds - hours*60*60 - minutes*60;
disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])

%% Stuff for figures
wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
wavemovie2(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);

f = figure();
subplot(2,1,1)
plot(xx, Dwave1.*conj(Dwave1))
%axis([-10e-10,10e-10,0,0.01])
title('ground')
subplot(2,1,2)
plot(xx, Dwave2.*conj(Dwave2))
%axis([-10e-10,10e-10,0,0.01])
title('excited')
wavemovie(1) = getframe(f);
close all
%{
f = figure();
subplot(2,1,1)
plot(xx, Awave1.*conj(Awave1))
title('ground')
subplot(2,1,2)
plot(xx, Awave2.*conj(Awave2))
title('excited')
wavemovie2(1) = getframe(f);
close all
%}
%{
g = figure();
plotThatStuff2(Dwave1, Dwave2, XX, YY, zlimit2, LL)
wavemovie2(1) = getframe(g);
close all
%}
%% Propagate

time = 0;
conservation = zeros(tNN,1);

tic

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,25) == 0
        
       
        f = figure();
        subplot(2,1,1)
        plot(xx, Dwave1.*conj(Dwave1))
        %axis([-10,10,0,0.01])
        title('ground')
        subplot(2,1,2)
        plot(xx, Dwave2.*conj(Dwave2))
        %axis([-10,10,0,0.01])
        title('excited')
        wavemovie(ll/25+1) = getframe(f);
        close all
        %{
        f = figure();
        subplot(2,1,1)
        plot(xx, Awave1.*conj(Awave1))
        title('ground')
        subplot(2,1,2)
        plot(xx, Awave2.*conj(Awave2))
        title('excited')
        wavemovie2(ll/25+1) = getframe(f);
        close all
        %}
        time = toc + time;
        tic
        averagetime = time/ll;
        
        seconds = (tNN-ll)*averagetime;
        hours = floor(seconds/60/60);
        minutes = floor((seconds - hours*60*60)/60);
        seconds = seconds - hours*60*60 - minutes*60;
        disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])
        
    end
    
    %% Hamiltonians
    
    W12 = -nu*xx.*sin(omegaLight*dt*ll);

    [V1, V2, A, Ainv] = makeW1D(W11, W12, W22);

    %V1 = V1 + exp(-(XX.^2+YY.^2))*1i;
    %V2 = V2 + exp(-(XX.^2+YY.^2))*1i;
    
    U1 = exp(-1i*V1/hbar*dt);
    U2 = exp(-1i*V2/hbar*dt);

    
    %% Propagation
    Awave1 = U1.*Awave1;
    Awave2 = U2.*Awave2;

    Dwave1 = (Ainv(:,1,1))'.*Awave1 + (Ainv(:,1,2))'.*Awave2;
    Dwave2 = (Ainv(:,2,1))'.*Awave1 + (Ainv(:,2,2))'.*Awave2;

    
    Dwave1 = ifft(KEUU.*fft(Dwave1));
    Dwave2 = ifft(KEUU.*fft(Dwave2));
    
    Awave1 = (A(:,1,1))'.*Dwave1 + (A(:,1,2))'.*Dwave2;
    Awave2 = (A(:,2,1))'.*Dwave1 + (A(:,2,2))'.*Dwave2;

    
    %Awave1 = U1.*Awave1;
    %Awave2 = U2.*Awave2;
    
    %Awave1 = Awave1/sqrt(sum(sum(Awave1.*conj(Awave1)))+sum(sum(Awave2.*conj(Awave2))));
    %Awave2 = Awave2/sqrt(sum(sum(Awave1.*conj(Awave1)))+sum(sum(Awave2.*conj(Awave2))));
    
    conservation(ll) = sum(sum(Dwave2.*conj(Dwave2)))+sum(sum(Dwave1.*conj(Dwave1)));
    
end
%}
%% Done :)

t = datetime;

%save('C:\Chelsea\TDSE_Output\A0eV.mat')

save(['C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])

%save(['C:\Chelsea\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])

x=linspace(0,2001,200);
y=sin(x);
sound(y)
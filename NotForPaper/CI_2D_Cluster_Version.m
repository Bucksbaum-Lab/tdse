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
parpool(24)

%% Discritize the wavefunction

LL = 50;
NN = 2^10;
xx = linspace(-LL, LL-2*LL/NN, NN); %Angstrom
dt = 1/100; %femtoseconds
totalTime = 10; %femtoseconds

dx = LL/NN;
tNN = floor(totalTime/dt);

kk = linspace(-1/dx, 1/dx, NN);

kk = fftshift(kk);

[XX,YY] = meshgrid(xx,xx);
[k1,k2] = meshgrid(kk,kk);

%% Constants

kappa = 1; %eV/Angstrom
lambda = 1; %eV/Angstrom
mm = 1; %amu
hbar = 1;
lighteV = 1; %eV
nu = 0; %eV/Angstrom

omegaLight = lighteV/4.13567;

%% Hamiltonians

tic
W11 = kappa*XX;
W12 = lambda*YY + nu*XX.*sin(omegaLight*0);
W22 = -kappa*XX;

[V1, V2, A, Ainv] = makeW(W11, W12, W22);

V1 = V1 + exp(-(XX.^2+YY.^2)/100000)*1i;
V2 = V2 + exp(-(XX.^2+YY.^2)/100000)*1i;

KEUU = exp(-1i*pi^2*hbar/8/mm*dt*(k1.^2+k2.^2)*0.0398);

U1 = exp(-1i*V1/hbar*dt/2*1.512);
U2 = exp(-1i*V2/hbar*dt/2*1.512);

%% Initial wavefunctions

Dwave1 = exp(-mm/2/hbar*((XX).^2+(YY).^2)/.1);
Dwave1 = Dwave1/sqrt(sum(sum(Dwave1.*Dwave1)));
Dwave2 = zeros(size(Dwave1));

Awave1 = A(:,:,1,1).*Dwave1 + A(:,:,1,2).*Dwave2;
Awave2 = A(:,:,2,1).*Dwave1 + A(:,:,2,2).*Dwave2;

zlimit = max(max(abs(Dwave1))).*1.1;
zlimit2 = max([sum(abs(Dwave1),1),(sum(abs(Dwave1),2))'])*1.1;

seconds = tNN*toc;
hours = floor(seconds/60/60);
minutes = floor((seconds - hours*60*60)/60);
seconds = seconds - hours*60*60 - minutes*60;
disp([num2str(hours) ' hours ' num2str(minutes) ' minutes ' num2str(seconds) ' seconds left'])

%% Stuff for figures
%{
wavemovie(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);
wavemovie2(floor(tNN/25)+1) = struct('cdata',[],'colormap',[]);

f = figure();
plotThatStuff(Dwave1, Dwave2, XX, YY, zlimit, LL)
wavemovie(1) = getframe(f);
close all

g = figure();
plotThatStuff2(Dwave1, Dwave2, XX, YY, zlimit2, LL)
wavemovie2(1) = getframe(g);
close all
%}

output1{1} = Dwave1;
output2{1} = Dwave2;
%% Propagate

time = 0;
conservation = zeros(tNN,1);

tic

for ll = 1:tNN
 
    %% Plotting stuff
    if mod(ll,25) == 0
        
        %{
        f = figure();

        plotThatStuff(Dwave1, Dwave2, XX, YY, zlimit, LL)
        
        wavemovie(ll/25+1) = getframe(f);
        close all

        g = figure();

        plotThatStuff2(Dwave1, Dwave2, XX, YY, zlimit2, LL)
        
        wavemovie2(ll/25+1) = getframe(g);
        close all
        %}
        
        output1{ll/25+1} = Dwave1;
        output2{ll/25+1} = Dwave2;
        
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
    
    W12 = lambda*YY + nu*XX.*sin(omegaLight*dt*ll);

    [V1, V2, A, Ainv] = makeW(W11, W12, W22);

    V1 = V1 + exp(-(XX.^2+YY.^2)/100000)*1i;
    V2 = V2 + exp(-(XX.^2+YY.^2)/100000)*1i;
    
    U1 = exp(-1i*V1/hbar*dt/2);
    U2 = exp(-1i*V2/hbar*dt/2);

    
    %% Propagation
    Awave1 = U1.*Awave1;
    Awave2 = U2.*Awave2;

    Dwave1 = A(:,:,1,1).*Awave1 + A(:,:,1,2).*Awave2;
    Dwave2 = A(:,:,2,1).*Awave1 + A(:,:,2,2).*Awave2;
    
    Dwave1 = ifft2(KEUU.*fft2(Dwave1));
    Dwave2 = ifft2(KEUU.*fft2(Dwave2));
    
    Awave1 = Ainv(:,:,1,1).*Dwave1 + Ainv(:,:,1,2).*Dwave2;
    Awave2 = Ainv(:,:,2,1).*Dwave1 + Ainv(:,:,2,2).*Dwave2;
    
    Awave1 = U1.*Awave1;
    Awave2 = U2.*Awave2;
    
    Awave1 = Awave1/sqrt(sum(sum(Awave1.*conj(Awave1)))+sum(sum(Awave2.*conj(Awave2))));
    Awave2 = Awave2/sqrt(sum(sum(Awave1.*conj(Awave1)))+sum(sum(Awave2.*conj(Awave2))));
    
    conservation(ll) = sum(sum(Dwave2.*conj(Dwave2)))+sum(sum(Dwave1.*conj(Dwave1)));
    
end

%% Done :)
%x=linspace(0,2001,200);
%y=sin(x);
%sound(y)

t = datetime;

save(['~/tdseOutput/Output' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])

%save('C:\Chelsea\TDSE_Output\A0eV.mat')

%save('C:\Users\chelsea\OneDrive\Documents\PhD\Acetylene\TDSE_Output\10eV.mat')

%save(['C:\Chelsea\TDSE_Output\' num2str(day(t)) '-' num2str(month(t)) '-' num2str(year(t)) '_' num2str(hour(t)) '-' num2str(minute(t)) '.mat'])
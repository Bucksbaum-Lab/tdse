close all
clear all

% Metric constants

% Hartree atomic units
hbar=1;
%numbers
%angstroms
L = 10;
fs = L/2^9;
m = 1;
N = L/fs;
nmax = 100;
ts = 1;

%x space, momentum space, and time
x = -L/2:fs:L/2-fs;
[x1,x2] = meshgrid(x,x);
k=(-N/2:N/2-1)/L;
%k = fftshift(k);
[k1,k2]=meshgrid(k,k);
%k1 = fftshift(fftshift(k1).').';
%k2 = fftshift(fftshift(k2).').';
tgrid = linspace(0,ts*nmax,nmax+1);

%build hamiltonian
KE = exp(hbar^2/(2*m)*(-(2*pi*k1).^2/(1i*hbar)*ts)+hbar^2/(2*m)*(-(2*pi*k2).^2/(1i*hbar)*ts));
%KE = ones(size(k1));

%Potential numbers
delta = 1;
lambda = 1;
kappa1 = 0;
kappa2 = 0;
absorption_strength = 1;
boundry_size = 2*L;

%electric field params
mu = 1;
E0 = 1;
omega = 1;

%W11 = -delta*x1;
%W22 = delta*x1;
%W12 = @(t)lambda*x2;
%+mu*E0*x1./sqrt(x1.^2+x2.^2+.01)*sin(omega*t);

%VG = @(t)exp(-(.5*(W11+W22)+sqrt((.5*(-W22+W11)).^2+W12(t).^2))/(1i*hbar)*ts);
%VE = @(t)exp(-(.5*(W11+W22)-sqrt((.5*(-W22+W11)).^2+W12(t).^2))/(1i*hbar)*ts);

%VG = -1./(sqrt((x1.^2+x2.^2)/1+1));
%VB = -1i*absorption_strength*(sin(pi*x1/(L))).^2-1i*absorption_strength.*(sin(pi*x2/(L))).^2;

%VG = exp(-x1*0/(1i*hbar)*ts);
%VE = exp(-x1*0/(1i*hbar)*ts);

%VG = exp((VG)/(1i*hbar)*ts);

%build initial wavefunction
%psi0 = (exp(-(x1).^2/20-(x2).^2/20));
%psi0 = psi0/sqrt(sum(sum(sum(psi0.^2*2*L*fs))));
wavefunctionG = rand(size(x1))+rand(size(x1))*1i;
wavefunctionG = wavefunctionG/sqrt(sum(sum(wavefunctionG.*conj(wavefunctionG))));

%wavefunctionE = psi0;

Mag(nmax+1) = struct('cdata',[],'colormap',[]);
Angle(nmax+1) = struct('cdata',[],'colormap',[]);

figure()
%subplot(2,1,1)
pcolor(x1,x2,sqrt(wavefunctionG.*conj(wavefunctionG)))
xlabel('x1')
ylabel('x2')
shading flat
title('ground mag')
colorbar
axis('equal')
axis([-L/2,L/2,-L/2,L/2])
%{
subplot(2,1,2)
pcolor(x1,x2,sqrt(wavefunctionE.*conj(wavefunctionE)))
xlabel('x1')
ylabel('x2')
shading flat
title('excited mag')
colorbar
axis('equal')
axis([-L/2,L/2,-L/2,L/2])
%}
drawnow

Mag(1) = getframe(gcf);

figure()
%subplot(2,1,1)
pcolor(x1,x2,unwrap(angle(wavefunctionG)))
xlabel('x1')
ylabel('x2')
shading flat
title('ground phase')
colorbar
axis('equal')
axis([-L/2,L/2,-L/2,L/2])
%{
subplot(2,1,2)
pcolor(x1,x2,unwrap(angle(wavefunctionE)))
xlabel('x1')
ylabel('x2')
shading flat
title('excited phase')
colorbar
axis('equal')
axis([-L/2,L/2,-L/2,L/2])
%}
drawnow
Angle(1) = getframe(gcf);
%}
%Adiabatic to diabatic transform

%AA11 = cos(.5*atan(2*W12(0)./(W11-W22)));
AA11 = ones(size(wavefunctionG));
%AA12 = sin(.5*atan(2*W12(0)./(W11-W22)));
AA12 = 0;
%AA11(x1==0&x2==0) = cos(.5*atan(lambda/(delta)));
%AA12(x1==0&x2==0) = sin(.5*atan(lambda/(delta)));

DwaveG = AA11.*wavefunctionG;
%- AA12.*wavefunctionE;
%DwaveE = AA12.*wavefunctionG + AA11.*wavefunctionE;

tic
for nn = 1:nmax
    
    %close all
    
    %AA11 = cos(.5*atan(2*W12(tgrid(nn))./(W11-W22)));
    AA11 = ones(size(wavefunctionG));
    %AA12 = sin(.5*atan(2*W12(tgrid(nn))./(W11-W22)));
    AA12 = 0;
    %AA11(x1==0&x2==0) = cos(.5*atan(lambda/(delta)));
    %AA12(x1==0&x2==0) = sin(.5*atan(lambda/(delta)));
    
    momentumG = fft2(DwaveG);
    
    %{
    figure()
    pcolor(k1,k2,abs(momentumG))
    xlabel('k1')
    ylabel('k2')
    shading flat
    title('abs mom')
    colorbar
    axis('equal')
    %axis([-L/2,L/2,-L/2,L/2])
    
    figure()
    pcolor(k1,k2,angle(momentumG))
    xlabel('k1')
    ylabel('k2')
    shading flat
    title('angle mom')
    colorbar
    axis('equal')
    
    figure()
    pcolor(k1,k2,abs(KE.*momentumG))
    xlabel('k1')
    ylabel('k2')
    shading flat
    title('abs ke mom')
    colorbar
    axis('equal')
    
    figure()
    pcolor(k1,k2,angle(KE.*momentumG))
    xlabel('k1')
    ylabel('k2')
    shading flat
    title('angle ke mom')
    colorbar
    axis('equal')
    %}
    
    DwaveG = ifft2(KE.*momentumG);
    %momentumE = (fft2((DwaveE)'));
    %DwaveE = (ifft2((KE.*momentumE)));
    
    AwaveGp = AA11.*DwaveG;
    %+ AA12.*DwaveE;
    %AwaveEp = -AA12.*DwaveG + AA11.*DwaveE;
    
    %wavefunctionG=VG(tgrid(nn)).*AwaveGp;
    %wavefunctionE=VE(tgrid(nn)).*AwaveEp;
    
    wavefunctionG=VG.*AwaveGp;
    %wavefunctionE=VE.*AwaveEp;

    wavefunctionG = wavefunctionG/sqrt(sum(sum(wavefunctionG.*conj(wavefunctionG))));
    
    DwaveG = AA11.*wavefunctionG;
    %- AA12.*wavefunctionE;
    %DwaveE = AA12.*wavefunctionG + AA11.*wavefunctionE;
    
    %DwaveG = DwaveG/sum(sum(sqrt(abs(DwaveG).^2)));
    %DwaveE = DwaveE/sum(sum(sqrt(abs(DwaveE).^2)));
    
    %{
    figure()
    %subplot(2,1,1)
    pcolor(k1,k2,sqrt(real(momentumG.*conj(momentumG))))
    xlabel('x1')
    ylabel('x2')
    shading flat
    title('ground mom')
    colorbar
    axis('equal')
    %axis([-L/2,L/2,-L/2,L/2])
    
    subplot(2,1,2)
    pcolor(k1,k2,sqrt(real(momentumE.*conj(momentumE))))
    xlabel('x1')
    ylabel('x2')
    shading flat
    title('excited mom')
    colorbar
    axis('equal')
    %axis([-L/2,L/2,-L/2,L/2])
    
    drawnow
    pause(2)
    %}
    figure()
    %subplot(2,1,1)
    pcolor(x1,x2,sqrt(real(wavefunctionG.*conj(wavefunctionG))))
    xlabel('x1')
    ylabel('x2')
    shading flat
    title('ground mag')
    colorbar
    axis('equal')
    axis([-L/2,L/2,-L/2,L/2])
    %{
    subplot(2,1,2)
    pcolor(x1,x2,sqrt(real(wavefunctionE.*conj(wavefunctionE))))
    xlabel('x1')
    ylabel('x2')
    shading flat
    title('excited mag')
    colorbar
    axis('equal')
    axis([-L/2,L/2,-L/2,L/2])
    %}
    drawnow
    Mag(nn+1) = getframe(gcf);
    
    figure()
    %subplot(2,1,1)
    pcolor(x1,x2,unwrap(angle(wavefunctionG)))
    xlabel('x1')
    ylabel('x2')
    shading flat
    title('ground phase')
    colorbar
    axis('equal')
    axis([-L/2,L/2,-L/2,L/2])
    %{
    subplot(2,1,2)
    pcolor(x1,x2,unwrap(angle(wavefunctionE)))
    xlabel('x1')
    ylabel('x2')
    shading flat
    title('excited phase')
    colorbar
    axis('equal')
    axis([-L/2,L/2,-L/2,L/2])
    %}
    drawnow
    
    Angle(nn+1) = getframe(gcf);
    %}
end
toc

[h, w, p] = size(wavemovie(1).cdata);

hf = figure;
axis off
movie(hf, wavemovieA,1,2)
hf = figure;
axis off
movie(hf, wavemovie2A,1,2)
%}
%{
hf = figure;
axis off
movie(hf, wavemovieD,1,1)
%}
hf = figure;
axis off
movie(hf, anglemovieA,1,1)
figure()
plot(0:dt:tNN*dt-dt, Excited, 0:dt:tNN*dt-dt, Ground, 0:dt:tNN*dt-dt, conservation)
axis([0 totalTime 0 1])
legend('Excited', 'Ground', 'Total')
figure()
plot(0:dt:tNN*dt-dt, Dipolet)

%{
kappa
lambda
nu
lighteV
mm/.0096
LL
log2(NN)
1/dt
%}
function plotThatStuff2(wavefunctionG, wavefunctionE, XX, YY, zlimit, LL)

set(gcf, 'Position', [10,50,900,900])

subplot(2,2,1)
plot(XX(1,:), sum(abs(wavefunctionE)))
shading flat
axis([-LL,LL,0,zlimit])
xlabel('x')
title('excited')

subplot(2,2,2)
plot(YY(:,1), sum(abs(wavefunctionE),2))
shading flat
axis([-LL,LL,0,zlimit])
xlabel('y')
title('excited')

subplot(2,2,3)
plot(XX(1,:), sum(abs(wavefunctionG)))
shading flat
axis([-LL,LL,0,zlimit])
xlabel('x')
title('ground')

subplot(2,2,4)
plot(YY(:,1), sum(abs(wavefunctionG),2))
shading flat
axis([-LL,LL,0,zlimit])
xlabel('y')
title('ground')

drawnow
function plotThatStuff(wavefunctionG, wavefunctionE, x1, x2, zlimit, LL)

set(gcf, 'Position', [10,50,900,900])

subplot(2,2,1)
surf(x1, x2, abs(wavefunctionE))
shading flat
axis([-LL,LL,-LL,LL,0,zlimit])
%caxis([0,zlimit])
xlabel('x')
ylabel('y')
title('excited')

bool = abs(wavefunctionE) > 0.001;
Angle = angle(wavefunctionE);
Angle = bool.*Angle;

subplot(2,2,2)
surf(x1, x2, Angle)
shading flat
axis([-LL,LL,-LL,LL,0,zlimit])
%caxis([0,zlimit])
xlabel('x')
ylabel('y')
title('excited')

subplot(2,2,3)
surf(x1, x2, abs(wavefunctionG))
shading flat
axis([-LL,LL,-LL,LL,0,zlimit])
%caxis([0,zlimit])
xlabel('x')
ylabel('y')
title('ground')

bool = abs(wavefunctionG) > 0.001;
Angle = angle(wavefunctionG);
Angle = bool.*Angle;

subplot(2,2,4)
surf(x1, x2, Angle)
shading flat
axis([-LL,LL,-LL,LL,0,zlimit])
%caxis([0,zlimit])
xlabel('x')
ylabel('y')
title('ground')

drawnow
figure('units', 'inches', 'position', [.5 .5 7.5 6.5])
plot(linspace(0,15,length(GroundA)),GroundA,'linewidth', 3, 'color',[0,112,184]/255)
hold on
plot(linspace(0,15,length(GroundA)),ExcitedA, 'linewidth',3, 'color', 'k')
title('800 nm control')
xlabel('time (fs)')
ylabel('adiabatic state population')
set(gca, 'fontsize', 18)
legend('ground', 'excited', 'location', 'east')
dt = 1/10;
omegaLight = 1/9.5460;

tNN = 10000;

for ll = 1:tNN

x(ll) = cos(omegaLight*dt*ll);

end
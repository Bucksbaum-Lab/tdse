function A = GetA(alpha, nn)

A = zeros(nn,1);

A(1) = besselj(0,alpha);

for ii = 2:nn
    A(ii) = 2*besselj(ii-1,alpha);
end
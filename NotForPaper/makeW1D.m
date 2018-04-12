function [Vp, Vm, A, Ainv] = makeW1D(W11, W12, W22)

ll = length(W11);

W = zeros(ll,2,2);
A = W;
Ainv = W;
Vp = zeros(ll,1);
Vm = zeros(ll,1);

for nn = 1:ll
    
    %W(nn,mm,:,:) = [W11(nn,mm), W12(nn,mm); W12(nn,mm), W22(nn, mm)];
    [A(nn,:,:), WD, Ainv(nn,:,:)] = eig([W11(nn), W12(nn); W12(nn), W22(nn)]);
    %AINV(nn,mm,:,:) = inv(squeeze(A(nn,mm,:,:)));
    Vp(nn) = WD(1);
    Vm(nn) = WD(4);
                
end
        
Vp = Vp';
Vm = Vm';

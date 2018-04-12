function [Vp, Vm, A, Ainv] = makeW(W11, W12, W22)

%{
ll = length(W11);

W = zeros(ll,ll,2,2);
A = W;
Ainv = W;

Vp = 0.5*(W11+W22)+sqrt((0.5*(W11-W22)).^2+(W12).^2);
Vm = 0.5*(W11+W22)-sqrt((0.5*(W11-W22)).^2+(W12).^2);

beta = 0.5*atan(W12./(0.5*(W11-W22)));

for nn = 1:ll
    for mm = 1:ll
        
        size([cos(beta), sin(beta); -sin(beta), cos(beta)])
        size(A(nn,mm,:,:))
        
        A(nn,mm,:,:) = [cos(beta), sin(beta); -sin(beta), cos(beta)];
        Ainv(nn,mm,:,:) = [cos(beta), -sin(beta); sin(beta), cos(beta)];
                
    end
end
%}
ll = length(W11);

W = zeros(ll,ll,2,2);
A = W;
Ainv = W;
Vp = zeros(ll,ll);
Vm = zeros(ll,ll);

parfor nn = 1:ll
    for mm = 1:ll
        
        %W(nn,mm,:,:) = [W11(nn,mm), W12(nn,mm); W12(nn,mm), W22(nn, mm)];
        [A(nn,mm,:,:), WD, Ainv(nn,mm,:,:)] = eig([W11(nn,mm), W12(nn,mm); W12(nn,mm), W22(nn, mm)]);
        %AINV(nn,mm,:,:) = inv(squeeze(A(nn,mm,:,:)));
        Vp(nn,mm) = WD(1);
        Vm(nn,mm) = WD(4);
                
    end
end
        

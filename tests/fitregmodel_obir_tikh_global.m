function [pass,maxerr] = test(opt)

% Check that OBIR works with Tikhonov regularization using the fnnls solver

rng(1)
t1 = linspace(0,3,100);
t2 = linspace(0,4,200);
r = linspace(2,7,100);
P = dd_gauss2(r,[3.7,0.5,0.5,4.3,0.3]);
K1 = dipolarkernel(t1,r);
K2 = dipolarkernel(t2,r);

V1 = K1*P + whitegaussnoise(t1,0.05);
V2 = K2*P + whitegaussnoise(t2,0.05);

Pfit = fitregmodel({V1,V2},{K1,K2},r,'tikhonov','aic','obir',true);

% Pass: OBIR fits the distribution
pass = all(abs(Pfit - P)<4e-1);
maxerr = max(abs(Pfit - P));
 
if opt.Display
    plot(r,P,'k',r,Pfit,'r')
    legend('truth','OBIR')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    grid on, axis tight, box on
end

end
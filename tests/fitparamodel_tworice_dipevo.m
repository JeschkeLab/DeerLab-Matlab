function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a two-Rician model

t = linspace(0,5,300);
r = linspace(1,6,300);
parIn = [3 0.4 0.6 4.5 0.3 0.4];
P = dd_rice2(r,parIn);
K = dipolarkernel(t,r);
S = K*P;
par0 = [2 0.1 1 0.6 0.2 1];

[~,Pfit] = fitparamodel(S,@dd_rice2,r,K,par0,'multistart',10);

%Pass: distance distribution is well fitted
pass = all(abs(Pfit - P) < 1e-5);

maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
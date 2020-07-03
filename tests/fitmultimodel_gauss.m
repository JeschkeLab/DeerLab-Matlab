function [pass,maxerr] = test(opt)

% Test fitmultimodel() using Gaussian models

t = linspace(0,2,300);
r = linspace(2,6,300);
paramtrue = [4 0.2 0.4 4 1 0.4 3 0.4 0.2];
P = dd_gauss3(r,paramtrue);
K = dipolarkernel(t,r);
S = K*P;
[Pfit] = fitmultimodel(S,K,r,@dd_gauss,5,'aicc');

% Pass: distribution is well fitted
pass = all(abs(Pfit - P) < 4e-1);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
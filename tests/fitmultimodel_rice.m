function [pass,maxerr] = test(opt)

% Test fitmultimodel() using Rician models

t = linspace(0,2,100);
r = time2dist(t);
paramtrue = [4 0.2 0.4 4 1 0.4 3 0.4 0.4];
P = dd_rice3(r,paramtrue);
K = dipolarkernel(t,r);
S = K*P;
Pfit = fitmultimodel(S,K,r,@dd_rice,5,'aicc','TolFun',1e-8);

% Pass: distribution is well fitted
pass = all(abs(Pfit - P) < 7e-2);

pass = all(pass);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
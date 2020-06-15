function [pass,maxerr] = test(opt)

% Check that fitmultimodel() works well with upper/lower bounds (2 parameters)

t = linspace(0,3,100);
r = time2dist(t);
InputParam = [4 0.1 0.35 5 0.4 0.65];
P = dd_gauss2(r,InputParam);

K = dipolarkernel(t,r);
S = K*P;

Pfit1 = fitmultimodel(S,K,r,@dd_gauss,3,'aic','Lower',[2 0.05],'Upper',[6 0.5]);
Pfit2 = fitmultimodel(S,K,r,@dd_gauss,3,'aic','Lower',[2 0.05]);
Pfit3 = fitmultimodel(S,K,r,@dd_gauss,3,'aic','Upper',[6 0.5]);

% Pass 1-3: all distributions are well fitted
pass(1) = any(abs(Pfit1 - P) < 7e-5);
pass(2) = any(abs(Pfit1 - P) < 7e-5);
pass(3) = any(abs(Pfit1 - P) < 7e-5);

pass = all(pass);

maxerr = max(abs(Pfit1 - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit1,r,Pfit2,r,Pfit3)
   legend('truth','fit 1','fit 2','fit 3')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
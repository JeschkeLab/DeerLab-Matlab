function [pass,maxerr] = test(opt)

% Test that confidence levels can be adjusted via option

rng(1)
t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.01);

par0 = [2 0.2];
[parFit,Pfit,cistruct] = fitparamodel(S,@dd_gauss,r,K,par0);

parCI1 = cistruct.ci(95);
parCI2 = cistruct.ci(50);
parCI3 = cistruct.ci(0.25);

% Pass 1-2: confidence intervals behave as expected
pass(1) = all(all(abs(parFit - parCI1.') > abs(parFit - parCI2.')));
pass(2) = all(all(abs(parFit - parCI2.') > abs(parFit - parCI3.')));

pass = all(pass);

maxerr = NaN;

 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
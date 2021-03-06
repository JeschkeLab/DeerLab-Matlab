function [pass,maxerr] = test(opt)

% Check that fitsignal() works with a dipolar evolution function

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,150);
P = dd_gauss(r,[4.5 0.6]);

lam = 0.4;
Bmodel = @(t,lam) bg_exp(t,0.4,lam);
K = dipolarkernel(t,r,lam,Bmodel);
V = K*P + whitegaussnoise(t,0.01);

par0 = {[],0.5,0.5};
lb = {[],0.2,0.2};
ub = {[],0.5,0.5};
[Vfit,Pfit] = fitsignal(V,t,r,'P',@bg_exp,@ex_4pdeer,par0,lb,ub);

% Pass 1: signal is well fitted
pass(1) = all(abs(Vfit - V) < 3e-2);
% Pass 2: distribution is well fitted
pass(2) = all(abs(Pfit - P) < 4e-1);

pass = all(pass);

maxerr = max(abs(Pfit - P));

if opt.Display
   subplot(211)
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
   subplot(212)
   plot(t,V,'k.',t,Vfit,'r')
   legend('data','fit')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end
function [pass,maxerr] = test(opt)

% Check that fitsignal() passes options well

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,200);
P = dd_gauss(r,[4.5 0.6]);
V = dipolarsignal(t,r,P,'noiselevel',0.01);

% Pass options of fitregmodel
[Vfit1,Pfit1] = fitsignal(V,t,r,'P','none','none','RegOrder',2);
% Pass options of selregparam
[Vfit2,Pfit2] = fitsignal(V,t,r,'P','none','none','Search','golden');
% Pass options of fitregmodel and selregparam
[Vfit3,Pfit3] = fitsignal(V,t,r,'P','none','none','RegOrder',2,'Search','golden');

% Pass 1-3: signal is well fitted
pass(1) = all(abs(Vfit1 - V) < 3e-2);
pass(2) = all(abs(Vfit2 - V) < 3e-2);
pass(3) = all(abs(Vfit3 - V) < 3e-2);
% Pass 4-6: distribution is well fitted
pass(4) = all(abs(Pfit1 - P) < 4e-1);
pass(5) = all(abs(Pfit2 - P) < 4e-1);
pass(6) = all(abs(Pfit3 - P) < 4e-1);

pass = all(pass);

maxerr = max(abs(Pfit1 - P));

if opt.Display
   subplot(211)
   plot(r,P,'k',r,Pfit1,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
   subplot(212)
   plot(t,V,'k.',t,Vfit1,'r')
   legend('data','fit')
   xlabel('t [\mus]')
   ylabel('V(t)')
   grid on, axis tight, box on
end

end
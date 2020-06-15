function [pass,maxerr] = test(opt)

% Test the cosntruction of a mixed model of Gaussians

t = linspace(0,3,200);
r = linspace(2,6,100);
parIn1 = [3 0.5];
P1 = dd_gauss(r,parIn1);
parIn2 = [4 0.5];
P2 = dd_gauss(r,parIn2);
P = 0.7*P2 + 0.3*P1;

mixedModel = mixmodels(@dd_gauss,@dd_gauss);
parInMix = [parIn1 0.3 parIn2 0.7];
Pmix = mixedModel(r,parInMix);

% Pass: the models have been mixed properly
pass = all(abs(Pmix - P) < 1e-8);

maxerr = max(abs(Pmix - P));
 
if opt.Display
   plot(r,P,'k',r,Pmix)
   legend('truth','mix')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end
end
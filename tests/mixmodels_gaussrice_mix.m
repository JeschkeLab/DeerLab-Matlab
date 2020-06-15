function [pass,maxerr] = test(opt)

% Test the cosntruction of a mixed model of different basis functions

t = linspace(0,3,200);
r = linspace(2,6,100);
parIn1 = [2.5 0.3];
P1 = dd_gauss(r,parIn1);
parIn2 = [3.5 0.3];
P2 = dd_gauss(r,parIn2);
parIn3 = [4.5 0.3];
P3 = dd_rice(r,parIn3);
P = 0.4*P2 + 0.3*P1 + 0.3*P3;

mixedModel = mixmodels(@dd_gauss,@dd_gauss,@dd_rice);
parInMix = [parIn1 0.3 parIn2 0.4 parIn3 0.3];
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
function [pass,maxerr] = test(opt)

% Check that fitparamodel rescaling does not change the results

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.5]);
B = bg_exp(t,0.3);
K = dipolarkernel(t,r,0.3,B);

scale = 1e9;
V = K*P + whitegaussnoise(t,0.005);
par0 = [0.5 0.5 2 0.3];

lower = [0 0 1 0.1];
upper = [1 1 20 5];
mymodel = @(t,param)dipolarkernel(t,r,param(1),bg_exp(t,param(2)))*dd_gauss(r,param(3:4));

[parfit1,Vfit1] = fitparamodel(V*scale,mymodel,t,par0,lower,upper,'Rescale',true,'MultiStart',5);
[parfit2,Vfit2] = fitparamodel(V,mymodel,t,par0,lower,upper,'Rescale',false,'MultiStart',5);

Pfit1 = dd_gauss(r,parfit1(3:4));
Pfit2 = dd_gauss(r,parfit2(3:4));

%Pass 1-2: distance distributions are well fitted
pass(1) = all(abs(Pfit1 - P) < 2e-1);
pass(2) = all(abs(Pfit2 - Pfit1) < 1e-2);

%Pass 3-4: dipolar signals are well fitted
pass(3) = all(abs(Vfit1 - V*scale) < scale*2e-2);
pass(4) = all(abs(Vfit2 - V) < 2e-2);

pass = all(pass);

maxerr = max(abs(Pfit1 - P));

if opt.Display
   plot(r,P,'k',r,Pfit1,'r',r,Pfit2,'b')
   legend('truth','rescaled','normalized')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
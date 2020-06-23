function [pass,maxerr] = test(opt)

% Test all options of backgroundstart()

t = linspace(0,3.2,200);
r = time2dist(t);
S = dipolarkernel(t,r)*dd_gauss(r,[3,0.5]);
B = bg_exp(t,0.5);
lam = 0.5;
F = (1 - lam) + lam*S;
V = F.*B;

%Check whether control statements on too smal start/end values work
backgroundstart(V,t,@bg_exp,'SearchStart',0.01,'SearchEnd',0.1);
%Check other options
backgroundstart(V,t,@bg_exp,'EndCutoff',100);

%Control using the defaults
t0 = backgroundstart(V,t,@bg_exp,'SearchStart',0.1*max(t),'SearchEnd',0.6*max(t));
[Bfit] = fitbackground(V,t,@bg_exp,t0);

% Pass: the background fitted using the optimized start fits well
pass = abs(B - Bfit) < 1e-2;
maxerr = max(abs(B - Bfit));
 
%Plot if requested
if opt.Display
    plot(t,B,t,Bfit)
    xlabel('t [\mus]')
    ylabel('B(t) ')
    legend('Truth','Fit')
    axis tight, grid on, box on
end

end
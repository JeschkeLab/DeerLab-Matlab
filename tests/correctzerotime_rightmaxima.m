function [pass,maxerr] = test(opt)

% Check that zero-time correction works when maximum is at later times

t_true = linspace(-5,1,400);
r = time2dist(t_true);
P = dd_gauss(r,[4,0.2]);
V = dipolarkernel(t_true,r)*P;
tshift = 1.2343;

t = t_true + tshift;
[t_corr,czt] = correctzerotime(V,t);

% Pass 1: corrected time-axis is equal to original
pass(1) = all(abs(t_corr - t_true.') < 1e-10);
% Pass 2: zero-time is returned properly
pass(2) = abs(czt' - tshift) < 1e-10;

pass = all(pass);
 
maxerr = max(abs(t_corr - t.'));

if opt.Display
   plot(t,V,t_corr,V)
   xlabel('t [\mus]')
   ylabel('V')
   grid on, axis tight, box on
   legend('Input','Output')
end

end
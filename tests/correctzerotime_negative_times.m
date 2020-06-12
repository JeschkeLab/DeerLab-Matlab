function [pass,maxerr] = test(opt)

% Check that zero-time correction works with negative times

t_true = linspace(-2,5,400);
r = linspace(1,7,400);
P = dd_gauss(r,[4,0.2]);
V = dipolarkernel(t_true,r)*P;

tshift = 1.21343;
t = t_true + tshift;

[ct,czt] = correctzerotime(V,t);

% Pass 1: corrected time-axis is equal to original
pass(1) = all(abs(ct - t_true.') < 1e-10);
% Pass 2: zero-time is returned properly
pass(2) = abs(czt - tshift) < 1e-10;

maxerr = max(abs(ct - t.'));

if opt.Display
    plot(t,V,ct,V)
    xlabel('t [\mus]')
    ylabel('V(t)')
    grid on, axis tight, box on
    legend('Input','Output')
end

end
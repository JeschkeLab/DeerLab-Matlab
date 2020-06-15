function [pass,maxerr] = test(opt)

% Check dimensionality, non-negativity, and values at boundaries of model

info = dd_circle();

r = linspace(0,30,301);
par0 = [info.Start];

lower = [info.Lower];
upper = [info.Upper];

P1 = dd_circle(r,par0);
P2 = dd_circle(r.',par0);
P3 = dd_circle(r,lower);
P4 = dd_circle(r,upper);

% Pass 1: dimensionality is correct
pass(1) = isequal(P1,P2);
% Pass 2: non-negativity of default values fulfilled
pass(2) = all(P1 >= 0);
% Pass 3: non-negativity of default boundaries fulfilled
pass(3) = all(P1 >= 0) & all(P2 >= 0);
% Pass 4: there are no NaN values
pass(4) = all(~isnan(P1)) & all(~isnan(P2)) & all(~isnan(P3)) & all(~isnan(P4));
% Pass 5: integral is 1
dr = mean(diff(r));
pass(5) = abs(sum(P1)*dr-1)<1e-10;
% Pass 6: check range
r0 = par0(1);
R = par0(2);
pass(6) = all(P1(r<r0-R)==0) && all(P1(r>r0+R)==0);

% Pass 7: compatible with non-uniform distance vectors
rnus = sqrt(linspace(1.5,6^2,800));
P5 = dd_circle(rnus,par0);
pass(7) = round(trapz(rnus,P5),4) == 1;

pass = all(pass);

maxerr = NaN;

end
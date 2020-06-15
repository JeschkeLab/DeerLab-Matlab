function [pass,maxerr] = test(opt)

% Check dimensionality, non-negativity, and values at boundaries of model

info = dd_rice3();

r = linspace(0,50,500);
par0 = [info.Start];

lower = [info.Lower];
upper = [info.Upper];

P1 = dd_rice3(r,par0);
P2 = dd_rice3(r.',par0);
P3 = dd_rice3(r,lower);
P4 = dd_rice3(r,upper);

% Pass 1: dimensionality is correct
pass(1) = isequal(P1,P2);
% Pass 2: non-negativity of default values fulfilled
pass(2) = all(P1 >= 0);
% Pass 3: non-negativity of default boundaries fulfilled
pass(3) = all(P1 >= 0) & all(P2 >= 0);
% Pass 4: there are no NaN values
pass(4) = all(~isnan(P1)) & all(~isnan(P2)) & all(~isnan(P3)) & all(~isnan(P4));
% Pass 5: compatible with non-uniform distance vectors
rnus = sqrt(linspace(1.5,6^2,800));
P5 = dd_rice3(rnus,par0);
pass(5) = round(trapz(rnus,P5),4) == 1;

pass = all(pass);
 
maxerr = NaN;

end
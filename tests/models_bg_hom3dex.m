function [pass,maxerr] = test(opt)

% Basic functionality test of background model

model = @bg_hom3dex;
info = model();

t = linspace(-5,5,500);
par0 = [info.Start];

lower = [info.Lower];
upper = [info.Upper];

B1 = model(t,par0);
B2 = model(t.',par0);
B3 = model(t,lower);
B4 = model(t,upper);
 
t0 = 2.5;
B5 = model(t0,par0);
Bval0 = 0.882896642393692;

B6 = model(t,par0,1);
B7 = model(t,par0.*[2 1],0.5);

% Pass 1-2: dimensionality is correct
pass(1) = isequal(B1,B2);
pass(2) = iscolumn(B1) && iscolumn(B2);
% Pass 3: there are no NaN values
pass(3) = all(~isnan(B1)) & all(~isnan(B2)) & all(~isnan(B3)) & all(~isnan(B4));
% Pass 4: specific value is reproducible
pass(4) = abs(B5-Bval0) < 1e-6;
% Pass 5: pathway amplitude can be specified correctly
pass(5) = isequal(B6,B7);

pass = all(pass);

maxerr = abs(B5-Bval0);

end

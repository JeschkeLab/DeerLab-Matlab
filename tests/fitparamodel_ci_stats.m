function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a wormchain model
t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P;

par0 = [2 0.2];
[~,Pfit,cistruct] = fitparamodel(S,@dd_gauss,r,K,par0,'solver','lsqnonlin');

parmean = cistruct.mean;
parmedian = cistruct.median;
parp50 = cistruct.percentile(0.50);

% Pass 1-2: all statistical indicators should be equal
pass(1) = all(abs(parmean - parmedian)<1e-5);
pass(2) = all(abs(parmedian - parp50)<1e-5);

pass = all(pass);
maxerr = max(abs(parmean - parmedian));

end
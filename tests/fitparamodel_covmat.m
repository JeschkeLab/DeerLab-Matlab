function [pass,maxerr] = test(opt)

% Check that covariance matrix can be manually specified

t = linspace(0,8,300);
r = linspace(1,8,300);
parin  = [4.5,0.5];
P = dd_gauss(r,parin);

sig = 0.01;
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,sig);

covmat = sig^2*eye(numel(t));

[~,Pfit,paruq] = fitparamodel(S,@dd_gauss,r,K,'Covariance',covmat);

% Pass: it does not crash
pass = true;

maxerr = max(abs(Pfit - P));


end
function [pass,maxerr] = test(opt)

% Check that multiple covariance matrices can be manually specified

t = linspace(0,8,300);
r = linspace(1,8,300);
parin  = [4.5,0.5];
P = dd_gauss(r,parin);

sig1 = 0.01;
sig2 = 0.02;
K = dipolarkernel(t,r);

S1 = K*P + whitegaussnoise(t,sig1);
S2 = K*P + whitegaussnoise(t,sig2);

covmat1 = sig1^2*eye(numel(t));
covmat2 = sig2^2*eye(numel(t));

[~,Pfit,paruq] = fitparamodel({S1,S2},@dd_gauss,r,{K,K},'Covariance',{covmat1,covmat2});

% Pass: it does not crash
pass = true;

maxerr = max(abs(Pfit - P));


end
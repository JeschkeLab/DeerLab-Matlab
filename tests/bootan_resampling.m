function [pass,maxerr] = test(opt)

% Check that the different bootan resampling methods work

rng(1)

sig = 0.05;

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.8]);
K = dipolarkernel(t,r);
V = K*P;

parfit = fitparamodel(V,@dd_gauss,r,K);

Vfit = K*dd_gauss(r,parfit) + whitegaussnoise(t,sig);

[ci] = bootan(@bootfcn,V,Vfit,10,'resampling','residual');
[ci2] = bootan(@bootfcn,V,Vfit,10,'resampling','gaussian');

% Pass 1-2: both resampling methods yield similar results
pass(1) = all(abs(ci{1}.mean - ci2{1}.mean)<1e-1);
pass(2) = all(abs(ci{2}.mean - ci2{2}.mean)<1e-1);

pass = all(pass);

maxerr = max(abs(ci{2}.mean - ci2{2}.mean));

    function [parfit,Pfit] = bootfcn(V)
        parfit = fitparamodel(V,@dd_gauss,r,K);
        Pfit = dd_gauss(r,parfit);
    end


end
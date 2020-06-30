function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a wormchain model
t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.02);

par0 = [2 0.2];
[~,~,cistruct] = fitparamodel(S,@dd_gauss,r,K,par0,'solver','lsqnonlin');

Prmean = cistruct.pardist(1).pdf;
Pfwhm = cistruct.pardist(2).pdf;

% Distribution for <r>
xmean = cistruct.mean(1);
sig = cistruct.std(1);
x1 = linspace(xmean-4*sig,xmean+4*sig,500);
Prmean_ref = 1/sig/sqrt(2*pi)*exp(-((x1-xmean)/sig).^2/2);
Prmean_ref = Prmean_ref/sum(Prmean_ref);

% Distribution for FWHM
xmean = cistruct.mean(2);
sig = cistruct.std(2);
x2 = linspace(xmean-4*sig,xmean+4*sig,500);
Pfwhm_ref = 1/sig/sqrt(2*pi)*exp(-((x2-xmean)/sig).^2/2);
Pfwhm_ref = Pfwhm_ref/sum(Pfwhm_ref);

% Pass 1-2: all parameter distributions are returned correctly
pass(1) = all(abs(Prmean - Prmean_ref)<1e-5);
pass(2) = all(abs(Pfwhm - Pfwhm_ref)<1e-5);

pass = all(pass);

maxerr = max(abs(Prmean - Prmean_ref));

if opt.Display
    subplot(121)
   plot(x1,Prmean_ref,'k',x1,Prmean,'r')
   legend('truth','fit')
   xlabel('<r>')
   ylabel('PDF')
   grid on, axis tight, box on
   
   
    subplot(122)
   plot(x2,Pfwhm_ref,'k',x2,Pfwhm,'r')
   legend('truth','fit')
   xlabel('FWHM')
   ylabel('PDF')
   grid on, axis tight, box on
end

end
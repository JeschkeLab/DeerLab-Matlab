function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a wormchain model
t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.05);

par0 = [2 0.2];
[parFit,Pfit,cistruct] = fitparamodel(S,@dd_gauss,r,K,par0,'solver','lsqnonlin');

%Get 95% confidence intervals
parCI = cistruct.ci(95);
parFit = parFit(:);

%Propagate error to distributions
lb = zeros(numel(r),1);
Puq = cistruct.propagate(@(parfit)dd_gauss(r,parfit),lb,[]);
Puq = Puq.ci(95);

% Pass 1: confidence intervals behave as expected
pass(1) =  all(parFit < parCI(:,2)) & all(parFit > parCI(:,1));
% Pass 2: error is well propagated
pass(2) = all(Pfit < Puq(:,2)) & all(Pfit > Puq(:,1));

pass = all(pass);

maxerr = max(abs(Pfit - P));
 
neif opt.Display
    hold on
   plot(r,P,'k',r,Pfit,'r')
   fill([r fliplr(r)],[Puq(:,1); flipud(Puq(:,2))],'r','FaceAlpha',0.3,'LineStyle','none')
   hold off
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
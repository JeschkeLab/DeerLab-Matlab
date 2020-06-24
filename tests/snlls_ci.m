function [pass,maxerr] = test(opt)

% Test snlls to solve a nonlinear-constrained  + linear-regularized

% Prepare test data
r = linspace(0,8,200);
t = linspace(0,4,200);
lambda = 0.25;
K = dipolarkernel(t,r,lambda);
parin = [3.5 0.4 0.6 4.5 0.5 0.4];
P = dd_gauss2(r,parin);
V = K*P;

% Non-linear parameters
% nlpar = [lam];
nlpar0 = 0.2;
lb = 0;
ub = 1;
% Linear parameters: non-negativity
lbl = zeros(numel(r),1);
ubl = [];
% Separable LSQ fit
[lambdafit,linfit,puq] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl);

% Get model fit
Pfit = linfit(:);
Vfit = Kmodel(lambdafit,t,r)*Pfit;

%Request confdence intervals
Pci = puq.ci(95,'lin');
lambdaci = puq.ci(95,'nonlin');

% Pass 1: modulation depth is well fitted
pass(1) = all(abs(lambda - lambdafit)<1e-1);
% Pass 2: the distribution is well fitted
pass(2) = all(abs(P - Pfit)<1e-2);
% Pass 3: all fits are within the CIs
pass(3) = lambdafit>lambdaci(1) & lambdafit<lambdaci(2);  
pass(4) = all(Pfit>=Pci(:,1)) & all(Pfit<=Pci(:,2));  

pass = all(pass);

maxerr = max(abs(P - Pfit));

if opt.Display
    subplot(211)
    plot(t,V,'k.',t,Vfit);
    subplot(212)
    cla,hold on
    plot(r,P,'k',r,Pfit,'r')
    fill([r fliplr(r)],[Pci(:,1); flipud(Pci(:,2))],'r','FaceAlpha',0.2,'LineStyle','none')
    hold off
    axis tight, grid on
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    legend('truth','fit')
end

    function K = Kmodel(p,t,r)
        
        % Unpack parameter
        lam  = p;
        
        % Generate basic kernel
        K = dipolarkernel(t,r,lam);
        
    end

end


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
[nonlinfit,linfit] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl);

% Get model fit
Pfit = linfit(:);
Vfit = Kmodel(nonlinfit,t,r)*Pfit;

% Pass 1: modulation depth is well fitted
pass(1) = all(abs(lambda - nonlinfit)<1e-1);
% Pass 2: the distribution is well fitted
pass(2) = all(abs(P - Pfit)<1e-2);

pass = all(pass);

maxerr = max(abs(P - Pfit));

if opt.Display
    subplot(211)
    plot(t,V,'k.',t,Vfit);
    subplot(212)
    plot(r,P,'k',r,Pfit,'r')
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


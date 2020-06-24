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
[lambdafit,linfit,~,stats] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl);

% Get model fit
Pfit = linfit(:);
Vfit = Kmodel(lambdafit,t,r)*Pfit;

% Pass 1: goodness of fit statistics are computed correctly
pass(1) = isstruct(stats);

pass = all(pass);

maxerr = max(abs(P - Pfit));

if opt.Display
    subplot(211)
    plot(t,V,'k.',t,Vfit);
    subplot(212)
    cla,hold on
    plot(r,P,'k',r,Pfit,'r')
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


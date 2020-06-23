function [pass,maxerr] = test(opt)

% Test snlls to solve a nonlinear-constrained  + linear-constrained

% Prepare test data
r = linspace(0,8,200);
t = linspace(0,4,200);
K = dipolarkernel(t,r);
parin = [3.5 0.4 0.6 4.5 0.5 0.4];
P = dd_gauss2(r,parin);
V = K*P;

% Non-linear parameters
% nlpar = [r1 w1 r2 w2];
nlpar0 = [3.2 0.2 4.2 0.3];
lb = [1 0.1 1 0.1];
ub = [20 5 20 5];
% Linear parameters
lbl = [0 0];
ubl = [1 1];
% Separable LSQ fit
[nonlinfit,linfit] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl);

% Get model fit
Pfit = dd_gauss2(r,[nonlinfit([1 2]) linfit(1) nonlinfit([3,4]) linfit(2)]);

% Pass 1: all non-linear parameters are well fitted
pass(1) = all(abs(parin([1 2 4 5]) - nonlinfit)<1e-1);
% Pass 2: all linear parameters are well fitted
pass(2) = all(abs(parin([3 6]) - linfit)<1e-1);
% Pass 3: the distribution is well fitted
pass(3) = all(abs(parin([3 6]) - linfit)<1e-1);

pass = all(pass);

maxerr = max(abs(P - Pfit));

if opt.Display
    plot(r,P,'k',r,Pfit,'r')
    axis tight, grid on
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    legend('truth','fit')
end

    function K = Kmodel(p,t,r)
        
        % Unpack parameter
        r1  = p(1); r2  = p(3);
        w1  = p(2); w2  = p(4);
        
        % Generate basic kernel
        K = dipolarkernel(t,r);
        
        % Get Gauss basis functions
        P1 = dd_gauss(r,[r1,w1]);
        P2 = dd_gauss(r,[r2,w2]);
        
        % Combine all non-linear functions into one
        K = [K*P1 K*P2];
        
    end

end


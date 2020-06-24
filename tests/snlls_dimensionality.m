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
[nonlinfit1,linfit1] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl);
[nonlinfit2,linfit2] = snlls(V.',@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl);
[nonlinfit3,linfit3] = snlls(V,@(p)Kmodel(p,t,r),nlpar0.',lb,ub,lbl,ubl);
[nonlinfit4,linfit4] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb.',ub,lbl,ubl);
[nonlinfit5,linfit5] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub.',lbl,ubl);
[nonlinfit6,linfit6] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl.',ubl);
[nonlinfit7,linfit7] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl.');

% Pass 1-2: all outputs are equal
pass(1) = isequal(nonlinfit1,nonlinfit2,nonlinfit3,nonlinfit4,nonlinfit5,nonlinfit6,nonlinfit7);
pass(2) = isequal(linfit1,linfit2,linfit3,linfit4,linfit5,linfit6,linfit7);

pass = all(pass);

maxerr = NaN;

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


function [pass,maxerr] = test(opt)

% Test snlls input schemes

% Prepare test data
r = linspace(0,8,100);
t = linspace(0,4,100);
K = dipolarkernel(t,r);
parin = [3.5 0.4 0.6 4.5 0.5 0.4];
P = dd_gauss2(r,parin);
V = K*P;

% Non-linear parameters: [r1 w1 r2 w2];
nlpar0 = [1 1 2 2];
lb = [];
ub = [];
% Linear parameters
lbl = [];
ubl = [];

% Separable LSQ fit
[nonlinfit1,linfit1] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,ubl,'TolFun',1e-3);
[nonlinfit2,linfit2] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,lbl,'TolFun',1e-3);
[nonlinfit3,linfit3] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,ub,'TolFun',1e-3);
[nonlinfit4,linfit4] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,lb,'TolFun',1e-3);
[nonlinfit5,linfit5] = snlls(V,@(p)Kmodel(p,t,r),nlpar0,'TolFun',1e-3);

% Pass 2: all inputs schemes are equivalent
pass(1) = isequal(nonlinfit1,nonlinfit2,nonlinfit3,nonlinfit4,nonlinfit5);
pass(2) = isequal(linfit1,linfit2,linfit3,linfit4,linfit5);

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


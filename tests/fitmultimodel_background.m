function [pass,maxerr] = test(opt)

% Check that fitmultigauss() can fit the background and modulation depth

t = linspace(-0.3,4,100);
r = linspace(3,6,200);
InputParam = [4 0.2 0.5 4.3 0.3 0.4];
P = dd_gauss2(r,InputParam);
B = bg_exp(t,0.15);
V = dipolarsignal(t,r,P,0.25,B,'noiselevel',0.0);

ub = {[0.9 1],[6 1]};
lb = {[0.2 0.01],[1 0.1]};

function K = Kmodel(par)
    lam = par(1);
    k = par(2);
    B = bg_exp(t,k);
    K = dipolarkernel(t,r,lam,B);
end

[Pfit,parfit] = fitmultimodel(V,@Kmodel,r,@dd_gauss,2,'aicc',lb,ub,'multistart',10);
Kparfit = parfit{1};
Vfit = Kmodel(Kparfit)*Pfit;
% Pass: distribution is well fitted
pass = all(abs(Pfit - P) < 2e-1);

maxerr = max(abs(Pfit - P));

if opt.Display
    subplot(211)
    plot(t,V,'k.',t,Vfit)
    legend('data','fit')
    xlabel('t [\mus]')
    ylabel('V(t)')
    grid on, axis tight, box on
    
    subplot(212)
    plot(r,P,'k',r,Pfit,'r')
    legend('truth','fit')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    grid on, axis tight, box on
end


end
function [pass,maxerr] = test(opt)

% Test time-domain global fit of parametric models

rng(2)
dt = 0.008;
r = linspace(1,5,100);
P = dd_gauss2(r,[3,0.6,0.5,4,0.6,0.25]);

Ntime1 = 100;
t1 = linspace(0,dt*Ntime1,Ntime1);
K1 = dipolarkernel(t1,r);
S1 = K1*P + whitegaussnoise(Ntime1,0.01);

Ntime2 = 200;
t2 = linspace(0,dt*Ntime2,Ntime2);
K2 = dipolarkernel(t2,r);
S2 = K2*P + whitegaussnoise(Ntime2,0.02);

Ntime3 = 300;
t3 = linspace(0,dt*Ntime3,Ntime3);
K3 = dipolarkernel(t3,r);
S3 = K3*P + whitegaussnoise(Ntime3,0.02);

Ss = {S1,S2,S3};
info = dd_gauss2();

par0 = [2 0.1 0.5 5 0.1 0.5];
upper = [info.Upper];
lower = [info.Lower];

parglobal = fitparamodel(Ss,@myglobalmodel,{t1,t2,t3},par0,lower,upper);
Pglobal = dd_gauss2(r,parglobal);

% Pass: global fit is a better fit than the local one
pass = all(abs(Pglobal - P)<1e-2);
 
maxerr = max(abs(P - Pglobal));

if opt.Display
    subplot(121)
    hold on
    plot(t1,S1,'k')
    plot(t1,K1*Pglobal,'r')
    plot(t2,S2+1,'k')
    plot(t2,K2*Pglobal + 1,'r')
    plot(t3,S3+2,'k')
    plot(t3,K3*Pglobal + 2,'r')
    grid on, axis tight, box on
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    subplot(122)
    hold on
    plot(r,P,'k')
    plot(r,Pglobal,'r')
    legend('truth','global')
    grid on, axis tight, box on
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
end

    function V = myglobalmodel(t,par)
 
        K1 = dipolarkernel(t1,r);
        K2 = dipolarkernel(t2,r);
        K3 = dipolarkernel(t3,r);
        P = dd_gauss2(r,par);
        V{1} = K1*P;
        V{2} = K2*P;
        V{3} = K3*P;
    end

end
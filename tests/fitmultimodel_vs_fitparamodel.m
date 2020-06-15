function [pass,maxerr] = test(opt)

% Check that fitmultimodel() and fitparamodel() return the same answer

t = linspace(0,3,300);
r = linspace(2,6,300);
InputParam = [3 0.5 0.4 4 0.5 0.6];
P = dd_gauss2(r,InputParam);
K = dipolarkernel(t,r);
S = K*P;
par0 = [2 0.1 0.5 5 0.1 0.5];

[~,P_FP] = fitparamodel(S,@dd_gauss2,r,K,par0);
P_MG = fitmultimodel(S,K,r,@dd_gauss,6);

% Pass: fitparamodel and fitmultigauss find the same solution
pass = all(abs(P_FP - P_MG) < 1e-9);

maxerr = max(abs(P_FP - P_MG));
 
if opt.Display
   plot(r,P,'k',r,P_MG,r,P_FP)
   legend('truth','multigauss','fitparamodel')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end
function [pass,maxerr] = test(opt)

% Test if selectmethod can use start points specified by the user

t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss2(r,[3 0.3 0.5 4 0.3]);
K = dipolarkernel(t,r);
S = K*P;

models = {@dd_gauss,@dd_gauss2,@dd_gauss3};

par0{1} = [3 1];
par0{2} = [3 1 0.5 3 1];
par0{3} = [3 1 0.2 3 1 0.2 3 1];

[optimum,metric] = selectmodel(models,S,r,K,'aicc',par0,'Solver','lsqnonlin');

% Pass: the optimal model is found
pass = optimum==2;
 
maxerr = NaN;


if opt.Display
    plot(1:3,metric)
    set(gca,'xtick',[1 2 3],'xticklabel',{'One Rice','Two Rice','Three Rice'})
    xtickangle(45)
    ylabel('AICc')
    grid on, axis tight, box on
end

end


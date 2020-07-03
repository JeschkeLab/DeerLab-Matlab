%============================================================================
% DeerLab Example:
% Multi-Gauss fit of a 4-pulse DEER signal
%============================================================================

% This example showcases how to fit a simple 4-pulse DEER signal with
% background using a multi-Gauss model, i.e automatically optimizing the
% number of Gaussians in the model.

clear, clc, clf

% Generating a dataset
%-----------------------------------------------------------------------------

rng(1)
t = linspace(-0.25,4,300); % time axis, us
r = linspace(2.5,4.5,300); % distance axis, nm
param0 = [3 0.3 0.2 3.5 0.3 0.45 3.9 0.2 0.20]; % parameters for three-Gaussian model
P = dd_gauss3(r,param0); % ground truth distance distribution
lam = 0.3; % modulation depth
conc = 250; % spin concentration, uM
noiselvl = 0.005; % noise level

% Generate 4pDEER dipolar signal with noise
V = K4pdeer([lam,conc],t,r)*P + whitegaussnoise(t,noiselvl);

% Multi-Gauss fitting
%-----------------------------------------------------------------------------

% Parameter bounds:
%     lambda conc   rmean fwhm 
lb = {[  0   0.05],[  1   0.05]}; % lower bounds
ub = {[  1   1500],[  20   5  ]}; % upper bounds

% Prepare the kernel model
Kmodel = @(par) K4pdeer(par,t,r);
NGauss = 5; % maximum number of Gaussians

% Fit the kernel parameters with an optimized multi-Gauss distribution
[Pfit,param,Puq,paramuq,Nopt,metrics,Peval] = fitmultimodel(V,Kmodel,r,@dd_gauss,NGauss,'aic',lb,ub,'MultiStart',10);

% Extract the parameters
Kparfit = param{1};

% Get the time-domain fit
K = Kmodel(param{1});
Vfit = K*Pfit;

% Confidence intervals of the fitted distance distribution
Pci95 = Puq.ci(95); % 95%-confidence interval
Pci50 = Puq.ci(50); % 50%-confidence interval

% Akaike weights
%-----------------------------------------------------------------------------

% When comparing different parametric models is always a good idea to look
% at the Akaike weights for each model. They basically tell you the
% probability of a model being the most optimal choice.

% Compute the Akaike weights
dAIC = metrics - min(metrics);
Akaikeweights = 100*exp(-dAIC/2)/sum(exp(-dAIC/2));

% Plots
%-----------------------------------------------------------------------------
subplot(321), cla
hold on
plot(t,V,'k.','LineWidth',1)
plot(t,Vfit,'b','LineWidth',1.5)
plot(t,(1-Kparfit(1))*bg_hom3d(t,Kparfit(2),Kparfit(1)),'b--','LineWidth',1.5)
axis tight,grid on,box on
legend('data','Vfit','Bfit')
xlabel('t [\mus]'),ylabel('V(t)')
axis tight

subplot(322), cla
hold on
plot(r,P,'k','LineWidth',1.5)
plot(r,Pfit,'b','LineWidth',1.5)
fill([r fliplr(r)], [Pci50(:,1); flipud(Pci50(:,2))],'b','Linestyle','none','facealpha',0.45)
fill([r fliplr(r)], [Pci95(:,1); flipud(Pci95(:,2))],'b','Linestyle','none','facealpha',0.25)
axis tight,grid on,box on
legend('truth','optimal fit','95%-CI')
xlabel('r [nm]'), ylabel('P(r)')

subplot(323)
bar(metrics + abs(min(metrics)),'b','FaceAlpha',0.6)
axis tight,grid on,box on
ylabel('\Delta AIC')
xlabel('Number of Gaussians')

subplot(325)
bar(Akaikeweights,'b','FaceAlpha',0.6)
axis tight,grid on,box on
ylabel('Akaike Weight [%]')
xlabel('Number of Gaussians')

subplot(3,2,[4 6]), cla
hold on
for i=1:numel(Peval)
    plot(r,P + 2*i.','k',r,Peval{i} + 2*i.','b-','LineWidth',1.5)
end
axis tight,grid on,box on
set(gca,'ytick',2:2:2*NGauss,'yticklabel',1:NGauss)
xlabel('r [nm]')
ylabel('Number of Gaussians')
legend('truth','fit')

% Model function for a 4pDEER dipolar kernel 
%-----------------------------------------------------------------------------
function K = K4pdeer(par,t,r)

% Unpack parameters
lam = par(1);
conc = par(2);

% Simualte background
B = bg_hom3d(t,conc,lam);
% Generate dipolar kernel
K = dipolarkernel(t,r,lam,B);

end
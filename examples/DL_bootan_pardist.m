%============================================================================
% DeerLab Example:
% Bootstrapped distributions of fit parameters
%============================================================================

% This example shows how to generate probability density functions of
% values for fit parameters using bootstrapping, showcased for 5pDEER.

clear, clc, clf

% Generate data
%-----------------------------------------------------------------------------
rng(1)
t = -0.1:0.08:6.5;               % time axis, us
r = linspace(1.5,6,100);    % distance axis, ns
param0 = [3 0.3 0.2 3.5 0.3 0.65 3.8 0.2 0.15]; % parameters for three-Gaussian model
P = dd_gauss3(r,param0);         % model distance distribution
B = @(t,lam)bg_hom3d(t,300,lam); % background decay
exparam = [0.6 0.3 0.1 3.2];     % parameters for 5pDEER experiment
pathinfo = ex_5pdeer(exparam);   % pathways information

Vexp = dipolarsignal(t,r,P,pathinfo,B,'NoiseLevel',0.01);

% Analysis
%-----------------------------------------------------------------------------

% Run the fit once as usual, to check that the model fits the data
[exparfit,bgparfit,Vfit] = fit(Vexp,t,r);

% Bootstrapping with 100 samples
[bootuq] = bootan(@(V)fit(V,t,r),Vexp,Vfit,100,'Verbose',true);

% Extract the uncertainty quantification for the parameters
exparam_uq = bootuq{1};
bgparam_uq = bootuq{2};

% Extract distributions for the experiment parameters
Lam0_dist = exparam_uq.pardist(1);
lam1_dist = exparam_uq.pardist(2);
lam2_dist = exparam_uq.pardist(3);
T02_dist  = exparam_uq.pardist(4);

% Extract distributions for the background parameters
conc_dist = bgparam_uq.pardist(1);

% Plot
%-----------------------------------------------------------------------------
subplot(321)
fill(Lam0_dist.values,Lam0_dist.pdf,'b','FaceAlpha',0.4)
xline(exparfit(1),'k--','LineWidth',2)
axis tight,grid on
xlabel('\Lambda_0'),ylabel('PDF')
legend('Bootstrapped','Fit','Location','best')

subplot(322)
fill(lam1_dist.values,lam1_dist.pdf,'b','FaceAlpha',0.4)
xline(exparfit(2),'k--','LineWidth',2)
axis tight,grid on
xlabel('\lambda_1'),ylabel('PDF')

subplot(323)
fill(lam2_dist.values,lam2_dist.pdf,'b','FaceAlpha',0.4)
xline(exparfit(3),'k--','LineWidth',2)
axis tight,grid on
xlabel('\lambda_2'),ylabel('PDF')

subplot(324)
fill(T02_dist.values,T02_dist.pdf,'b','FaceAlpha',0.4)
xline(exparfit(4),'k--','LineWidth',2)
axis tight,grid on
xlabel('T_{0,2} [\mus]'),ylabel('PDF')

subplot(325)
fill(conc_dist.values,conc_dist.pdf,'b','FaceAlpha',0.4)
xline(bgparfit(1),'k--','LineWidth',2)
axis tight,grid on
xlabel('Spin conc. [\muM]'),ylabel('PDF')

function [exparam,bgparam,Vfit] = fit(V,t,r)

% Set boundaries for the fit parameters (see DL_fitting_5pdeer.m)
ex_lb   = [ 0    0    0   max(t)/2-1]; % lower bounds
ex_ub   = [10   10   10   max(t)/2+1]; % upper bounds
ex_par0 = [0.5  0.5  0.5  max(t)/2  ]; % start values
ub = {[],[],ex_ub};
lb = {[],[],ex_lb};
par0 = {[],[],ex_par0};

% Run the fit, since we are only interested in the parameters we'll ignore
% the rest (otherwise the Bfit,Pfit,etc. could be bootstrapped as well) 
% We need the Vfit to pass it to bootan as well, so we'll request that one too
[Vfit,~,~,parfit] = fitsignal(V,t,r,'P',@bg_hom3d,@ex_5pdeer,par0,lb,ub);

% Unpack the parameters, since bootan() requires the outputs to be arrays
% of numerical values, not structures
exparam = parfit.ex;
bgparam = parfit.bg;

end
%============================================================================
% DeerLab Example:
% Uncertainty propagation from parameter fits using covariance-based uncertainty
% quantificaion
%============================================================================

% This example explains how to propagate the uncertainty of the fitted
% parameters to the models which depend on them.

clc,clf,clear

% Generate data
%-----------------------------------------------------------------------------
rng(1)
t = linspace(-0.2,4,300);
r = linspace(2,5,400);
center = 3.5; % [nm] Rician center distance
width = 0.3; % [nm] Rician width
lam = 0.27; % Modulation depth
conc = 150; % [uM] Spin concentration
P = dd_rice(r,[center width]);
B = bg_hom3d(t,conc,lam);
V = dipolarsignal(t,r,P,lam,B,'noiselevel',0.03);

% Fit the data
%-----------------------------------------------------------------------------

% Pre-calculate the elemental dipolar kernel (for speed)
K0 = dipolarkernel(t,r);

% First we define the models for the different elements in our analysis
% (background, distribution and dipolar signal). For simplicity these
% models take the full parameter set
%
%          par = [lambda center width conc]; 
%
% and select the appropiate elements from the parameter set, i.e.
%
%   Pmodel = f1(center,width) -> par(2) & par(3)
%   Bmodel = f2(conc,lambda)  -> par(4) & par(1)
%   Vmodel = f3(par)          -> par(1) & par(2) & par(3) & par(4)
%
% By defining the models like this, we can spare then the indexing of the
% parameters each time we call one of these model and can pass the full
% parameter set directly.

Pmodel = @(par) dd_rice(r,par(2:3));
Bmodel = @(par) bg_hom3d(t,par(4),par(1));
Vmodel = @(par) (1 - par(1) + par(1)*K0*Pmodel(par)).*Bmodel(par);


% Next since we are dealing with a custom-defined model we need to specify
% the start values as well as boundaries of the parameter set:

% Parameters:[lam  center width  conc]
par0  =      [0.35  4.0   0.4    500 ]; % start values
lower =      [0.10  2.0   0.1    0.1 ]; % lower bounds
upper =      [0.50  7.0   0.5    1500]; % upper bounds

% Finally we can fit the data. It is important to request the parameter
% uncertainty quantification (paruq) output!
[fitpar,~,paruq] = fitparamodel(V,@(t,par)Vmodel(par),t,par0,lower,upper);

% Forward-calculate the models with the fitted parameters
Vfit = Vmodel(fitpar);
Pfit = Pmodel(fitpar);
Bfit = Bmodel(fitpar);
lamfit = fitpar(1);


% Uncertainty propagation
%-----------------------------------------------------------------------------

% In DeerLab, all uncertainty quantification structures contain a field
% .propagate() which is a function handle. This function has all the
% internal information on the covariance matrices required to propagate the
% uncertainty from the parameters to the models. 
%
% Thus, all we neeed to do is call .propagate and pass the model function
% which we want to propagate the uncertainty to. It is important that if
% the uncertainty quantification structure is defined for N-parameters (N=4
% in this case) the model function must accept all N parameters. Since we
% defined our model function to accept all N parameters already we do not
% need to worry about it.

% 1. Uncertainty of the dipolar signal fit
% =========================================
% This case is easy, we already have the model and it is unconstrained
Vuq = paruq.propagate(Vmodel); % Uncertainty quantification for Vfit
Vci95 = Vuq.ci(95); % 95%-confidence intervals for Vfit

% 2. Uncertainty of the distance distribution
% ===========================================
% In this case, the distribution has a non-negativity constraint which we
% can specify via the lb input. 
lb = zeros(numel(r),1); % Non-negativity constraint
Puq = paruq.propagate(Pmodel,lb); % Uncertainty quantification for Pfit
Pci95 = Puq.ci(95); % 95%-confidence intervals for Pfit

% 3. Uncertainty of the background
% ================================
% In this case, since we want to use this for plotting we need to evaluate
% the function (1-lambda)*Bfit instead of just Bfit in order to plot the\
% correct function.
Buq = paruq.propagate(@(p)(1-p(1))*Bmodel(p)); % Uncertainty quantification for (1-lam)Bfit
Bci95 = Buq.ci(95); % 95%-confidence intervals for (1-lam)Bfit

% Plots
%-----------------------------------------------------------------------------

% Time-domain
subplot(211),cla,hold on
plot(t,V,'k.',t,Vfit,'r',t,(1-lamfit)*Bfit,'b','LineWidth',1.5)
fill([t fliplr(t)],[Vci95(:,1); flipud(Vci95(:,2))],'r','FaceAlpha',0.3,'LineStyle','none')
fill([t fliplr(t)],[Bci95(:,1); flipud(Bci95(:,2))],'b','FaceAlpha',0.3,'LineStyle','none')
axis tight, grid on, box on
xlabel('t [\mus]'),ylabel('V(t)')
legend('data','Vfit','Bfit','Vfit 95%-CI','Bfit 95%-CI')

% Distance-domain
subplot(212),cla,hold on
plot(r,P,'k',r,Pfit,'r','LineWidth',1.5)
fill([r fliplr(r)],[Pci95(:,1); flipud(Pci95(:,2))],'r','FaceAlpha',0.3,'LineStyle','none')
xlabel('r [nm]'),ylabel('P(r) [nm^{-1}]')
axis tight, grid on, box on
legend('truth','Pfit','Pfit 95%-CI')

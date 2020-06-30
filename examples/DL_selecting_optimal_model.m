%============================================================================
% DeerLab Example:
% Selecting an optimal parametric model for fitting a dipolar signal
%============================================================================

% In this example we will show how to optimally select a parametric model for 
% a given dipolar signal.

clear,clc,clf

% Data Generation
%-----------------------------------------------------------------------------

% Let's start by constructing a simple dipolar signal with some noise arising 
% from a bimodal Gaussian distance distribution.

% Prepare the signal components
rng(1)
t = linspace(-0.3,6,300);
r = linspace(2,6,200);
P = dd_gauss2(r,[3.8 0.7 0.7 4.5 0.3 0.7]);

% Simulate the signal
V = dipolarsignal(t,r,P,'noiselevel',0.03);

% Plot
subplot(2,2,1)
plot(t,V,'k.')
axis tight, grid on
xlabel('t [\mus]'),ylabel('V(t)')


subplot(2,2,2)
plot(r,P,'k','Linewidth',1.5)
axis tight, grid on
xlabel('r [nm]'),ylabel('P(r) [nm^{-1}]')

% Selecting an optimal model
%-----------------------------------------------------------------------------

% Even though we know the ground truth, in this example we will cosider the 
% following set of potential parametric models: 
% - Unimodal Rician distribution
% - Bimodal Rician distribution
% - Trimodal Rician distribution
% - Unimodal Gaussian distribution
% - Bimodal Gaussian distribution
% - Trimodal Gaussian distribution
% - Mixed bimodal Gaussian/Rician distribution

% The first six models have built-in parametric models which we can use directly. 
% The last model we can construct from built-in models using the |mixmodels| function.

% Prepare the mixed model
dd_rice_gauss = mixmodels(@dd_rice,@dd_gauss);
 
% Prepare list of candidate parametric models
Models = {@dd_rice,@dd_rice2,@dd_rice3,...
          @dd_gauss,@dd_gauss2,@dd_gauss3,dd_rice_gauss};
     
% Prepare the dipolar kernel
K = dipolarkernel(t,r);

% Optimize the selection
[optIdx,AICc] = selectmodel(Models,V,r,K,'aicc');

dAICc = AICc - min(AICc);

% Plot the results
subplot(2,2,3)
cla,hold on
bar(dAICc,'k')
grid on, box on,axis tight
ylabel('\DeltaAIC')
tags = {'R_1','R_2','R_3','G_1','G_2','G_3','R/G'};
set(gca,'xtick',1:length(Models),'xticklabel',tags)
xlim([0 8]),ylim([0 0.7])

% Akaike Weights
%-----------------------------------------------------------------------------

% It is often more useful to look at these results from the perspective of
% Akaike weights, i.e. the probabilities of a model being the most optimal.

weights = 100*exp(-(dAICc/2))/sum(exp(-dAICc/2));

% Plot the results
subplot(2,2,4)
cla,hold on
bar(weights,'k')
grid on, box on,axis tight
ylabel('Akaike Weights [%]')
tags = {'R_1','R_2','R_3','G_1','G_2','G_3','R/G'};
set(gca,'xtick',1:length(Models),'xticklabel',tags)
xlim([0 8]),ylim([0 20])

% Typically there is not a single optimal model unless the noise level is very
% low. Usually several models have similar probabilities and should therefore be presented together. 



%============================================================================
% DeerLab Example:
% Extracting Gaussian constraints from a parameter-free distribution fit
%============================================================================

% While parameter-free distance distributions are the most robust way to
% analyze dipolar signals, many structural biology modelling programs
% accept only estimators such as mean distances or Gaussian constraints. 

% This example shows how to extract Gaussian constraints from a
% parameter-free fit of a dipolar signal and how to calculate the
% corresponding uncertainty. 

clf,clc,clear

% Generating a dataset
%-----------------------------------------------------------------------------

% For this example we will simulate a simple 4pDEER signal

% Parameters
t = linspace(0,5,250);
r = linspace(1,7,200);
P = dd_gauss3(r,[4.5 0.6 0.4 3 0.4 0.3 4 0.7 0.5]);
lambda = 0.3;
kappa = 0.4;

% Simulate the signal
Bmodel = @(t,lam) bg_hom3d(t,kappa,lam);
K = dipolarkernel(t,r,lambda,Bmodel);
V = K*P + whitegaussnoise(t,0.01);

% Fit the dipolar signal
%-----------------------------------------------------------------------------

% First, we need to fit the parameter-free distance distribution using fitsignal()
% We are only interested right now on the fitted distribution and the
% corresponding uncertainty quantification, so we will ignore the rest of
% the outputs.
[~,Pfit,~,~,fituq] = fitsignal(V,t,r,'P',@bg_exp,@ex_4pdeer,'Display',true);

% Extract Gaussian constraints from the fit
%-----------------------------------------------------------------------------

% Next, we will fit a multi-Gauss distribution to the fitted parameter-free
% distribution. We can do this by using the fitparamodel() function (in
% this example, fitting a two-Gauss model). 

% However, in order to get the correct uncertainty quantification, we need
% to specify the covariance matrix of the fitted distribution.
% fitparamodel() can then use that information to propagate the error in
% Pfit to the Gauss constraints that we then fit.

% Extract the uncertainty quantification of the fitted distribution...
Pfit_uq = fituq.Pfit;
% ...specifically its covariance matrix
Pfit_covmat = Pfit_uq.covmat;

% Fit a 2-Gauss model to the fitted parameter-free distribution:
%  - parfit: will contain the Gaussian constraints
%  - PGauss: the corresponding distribution
%  - paruq: the uncertainty quantification of our constraints
[parfit,PGauss,paruq] = fitparamodel(Pfit,@dd_gauss2,r,'Covariance',Pfit_covmat);

% Extract the 95%-confidence intervals...
par95 = paruq.ci(95);
% ... and print the results of the constraints 
fprintf('\nGaussian constraints:\n')
info = dd_gauss2();
for i=1:numel(parfit)
    fprintf('  parfit(%i) = %2.2f (%2.2f, %2.2f) %s\n',i,parfit(i),par95(i,1),par95(i,2),info(i).Parameter)
end

% Now propagate the error of the constraints on the model
lb = zeros(numel(r),1);% Non-negativity constraint
PGauss_uq = paruq.propagate(@(par)dd_gauss2(r,par),lb,[]);
PGauss95 = PGauss_uq.ci(95);

% Plot the fitted constraints model on top of the parameter-free case
hold on
plot(r,PGauss,'b','LineWidth',1.5)
fill([r fliplr(r)],[PGauss95(:,1); flipud(PGauss95(:,2))],'b','FaceAlpha',0.2,'LineStyle','none')
hold off
legend('Fit','95%-CI','50%-CI','2G-constraints','95%-CI')


%============================================================================
% DeerLab Example:
% Fitting a 5-pulse DEER signal with a parameter-free distribution
%============================================================================

% This example shows how to fit a 5-pulse DEER signal with a parameter-
% free distribution, a background, and all pathways parameters

clear, clc, clf

% Generate data
%-----------------------------------------------------------------------------
rng(1)
t = -0.1:0.02:6.5;               % time axis, us
r = linspace(1.5,6,numel(t));    % distance axis, ns
param0 = [3 0.3 0.2 3.5 0.3 0.65 3.8 0.2 0.15]; % parameters for three-Gaussian model
P = dd_gauss3(r,param0);         % model distance distribution
B = @(t,lam)bg_hom3d(t,300,lam); % background decay
exparam = [0.6 0.3 0.1 3.2];     % parameters for 5pDEER experiment
pathinfo = ex_5pdeer(exparam);   % pathways information

Vexp = dipolarsignal(t,r,P,pathinfo,B,'NoiseLevel',0.01);

% Run fit
%-----------------------------------------------------------------------------

% Now, 5pDEER data contain 3 additional parameters compared to 4pDEER (due
% to the additional dipolar pathway present in the signal). However, the
% refocusing time of the second dipolar pathway is very easy to constrain
% and strongly helps stabilizing the fit. 

% This pathway ususally refocuses at around t = max(t)/2, and usually can
% be even estimated from simple visual inspection of the signal. 
% Thus, we can strongly constraint this parameters while leaving the
% pathway amplitudes pretty unconstrained.

% Experiment parameters:
%          Lam0 lam1 lam2    T02
ex_lb   = [ 0    0    0   max(t)/2-1]; % lower bounds
ex_ub   = [100  100  100  max(t)/2+1]; % upper bounds
ex_par0 = [0.5  0.5  0.5  max(t)/2  ]; % start values

% In this case we only want to set the bounds for the experiment
% parameters, so we can leave the rest empty:
ub = {[],[],ex_ub};
lb = {[],[],ex_lb};
par0 = {[],[],ex_par0};

% Run the fit
[Vfit,Pfit,Bfit,parfit,parci] = ...
    fitsignal(Vexp,t,r,'P',@bg_hom3d,@ex_5pdeer,par0,lb,ub,'Display',true);

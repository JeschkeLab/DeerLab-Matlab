%
% SELECTMODEL Optimal parametric model selection
%
%   opt = SELECTMODEL({@model1,...,@modelN},V,r,K,{'aic',...})
%   Evaluates the fits of the parametric models (model1,...,modelN) to a
%   signal (V) according to the dipolar kernel (K) and distance axis (r).
%   The models must be passed as a cell array of function handles. Each fit
%   is then evaluated according to the model selection criterions
%   ('aic','aicc','bic','rmsd') specified in the last input argument M-point
%   cell array. Function returns a M-point array containing the optimal models
%   according to each selection method.
%
%   opt = SELECTMODEL({@model1,...,@modelN},V,t,{'aic',...})
%   Evaluates the fits of the time-domain parametric models by specifying
%   the time axis (t).
%
%   opt = SELECTMODEL({@model1,...,@modelN},{V1,V2,...},r,{K1,K2,...},{'aic',...})
%   Passing multiple signals/kernels enables distance-domain global fitting
%   of the parametric models to single distributions. 
%   The multiple signals are passed as a cell array of arrays of sizes N1,N2,...
%   and a cell array of kernel matrices with sizes N1xM,N2xM,... must be 
%   passed as well.
% 
%   opt = SELECTMODEL({@model1,...,@modelN},{V1,V2,...},{t1,t2,...},{'aic',...})
%   Similarly, time-domain global fitting can be used when passing time-domain
%   and the model time axes {t1,t2,...} of the corresponding signals.
%
%   opt = SELECTMODEL({@model1,...,@modelN},V,r,K,{'aic',...},{par1,...,parN})
%   The initial guess values for the parameters of each model can be passed
%   as a cell array {par1,...parN} of value vectors.
%
%   opt = SELECTMODEL({@model1,...,@modelN},V,r,K,{'aic',...},{par1,...,parN},{lb1,lb2,...,lbN},{ub1,ub2,...,ubN})
%   The upper and lower boundaries for the parameters of the individual
%   models can be specified as a cell array {lb1,lb2,...,lbN} and {ub1,ub2,...,ubN}
%   of value vectors 
%
%   [opt,f,params,paramcis,stats] = SELECTMODEL(__)
%   Returns a cell array the method selector functionals (f) for the
%   different methods and a cell array (params) with the fitted parameters
%   for each of the evaluated models as well as their confidence intervals (paramcis).
%   A cell array of structures (stats) containing the goodness-of-fit statistics 
%   of the different models are returned as well.
%
%   opt = SELECTMODEL(__,'Property',Value)
%   Additional (optional) arguments can be passed as name-value pairs.
%
%   'GlobalWeights' - M-element array of values used to weight the M input signals
%                     during global analysis
%
%   See "help fitparamodel" for a detailed list of other name-value pairs
%   accepted by the function.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [optima,functionals,fitparams,paramcis,stats] = selectmodel(models,Vs,ax,varargin)

if nargin<4
    error('At least four input arguments required.')
end

% Parse input schemes
%-------------------------------------------------------------------------------
% Prepare empty containers
Ks = [];
param0 = [];
methods = [];
Lower = [];
Upper = [];
allowedmethods = {'aic','aicc','rmsd','bic','all'};
%Check if kernel is passed
input = varargin{1};
isTimeDomain = true;
if ~ischar(input)
    if ~iscell(input), input = {input}; end
    if ~ischar(input{1})
        isTimeDomain = false;
        Ks = varargin{1};
        varargin(1) = [];
    end
end

% Parse the varargin cell array
optionstart = numel(varargin);
for i=1:numel(varargin)
    if ischar(varargin{i}) && all(~strcmp(varargin{i},allowedmethods))
        optionstart = i-1;
        break
    end
    switch i
        case 1
            methods = varargin{1};
        case 2           
            param0 = varargin{2};
        case 3
            Lower = varargin{3};
        case 4
            Upper = varargin{4};
    end
end
varargin(1:optionstart) = [];

if isempty(param0)
    param0 = {};
    param0(1:length(models)) = {[]};
end

% Input validation
%-------------------------------------------------------------------------------
if ~iscell(models)
    error('First input must be a cell array of model functions.');
end
if ~iscell(methods)
    methods = {methods};
end

% Parse the optional parameters in the varargin
[GlobalWeights] = parseoptions(varargin);

if ~isempty(Upper) && ~iscell(Upper)
    error('Upper property must be a cell array of upper bound vectors.')
end
if ~isempty(Lower) && ~iscell(Lower)
    error('Lower property must be a cell array of lower bound vectors.')
end
if length(Upper) > length(models) || length(Lower) > length(models)
    error('Lower/Upper bound cell array cannot exceed the number of models.')
end

%Parse the required inputs for global fitting
if ~iscell(Vs)
   Vs = {Vs}; 
end
if ~isTimeDomain && ~iscell(Ks)
   Ks = {Ks}; 
end
if isTimeDomain && ~iscell(ax)
   ax = {ax}; 
end
if isempty(GlobalWeights)
    GlobalWeights = globalweights(Vs);
else
    validateattributes(GlobalWeights,{'numeric'},{'nonempty','nonnegative'},mfilename,'GlobalWeights')
    GlobalWeights = GlobalWeights/sum(GlobalWeights);
end

% Set the bounds for the models
UpperBounds = cell(1,length(models));
LowerBounds = cell(1,length(models));
UpperBounds(1:length(Upper)) = Upper;
LowerBounds(1:length(Lower)) = Lower;


allowedMethodInputs = {'aic','aicc','bic','rmsd'};
for i = 1:length(methods)
    if strcmp(methods{i},'all')
        methods = allowedMethodInputs;
        break;
    end
    validateattributes(methods{i},{'char'},{'nonempty'})
    methods{i} = validatestring(methods{i},allowedMethodInputs);
end

% Run all parametric model fits and evaluate selection metrics
%-------------------------------------------------------------------------------
nMethods = length(methods);
nModels = length(models);
AICc = zeros(nModels,1);
BIC = zeros(nModels,1);
AIC = zeros(nModels,1);
RMSD = zeros(nModels,1);
fitparams = cell(nMethods,1);
paramcis = cell(nMethods,1);
stats = cell(nMethods,1);
for i = 1:nModels
    if isTimeDomain
        [parfit,fit,parci,~,stats{i}] = fitparamodel(Vs,models{i},ax,param0{i},LowerBounds{i},UpperBounds{i},varargin{:});
    else
        [parfit,fit,parci,~,stats{i}] = fitparamodel(Vs,models{i},ax,Ks,param0{i},LowerBounds{i},UpperBounds{i},varargin{:});
    end
    fitparams{i} = parfit;
    paramcis{i} = parci;
    
    if isTimeDomain && ~iscell(fit)
        fit = {fit};
    end
    
    nParams = numel(parfit);
    logprob = 0;
    for idx = 1:numel(Vs)
        V = Vs{idx};
        N = numel(Vs);
        if isTimeDomain
            SSR = sum((V(:) - fit{idx}).^2);
        else
            K = Ks{idx};
            SSR = sum((V(:) - K*fit).^2);
        end
        Q = nParams + 1;
        logprob = logprob + GlobalWeights(idx)*N*log(SSR/N);
        RMSD(i) = RMSD(i) + GlobalWeights(idx)*sqrt(1/N*SSR);
    end
    AIC(i) =  logprob + 2*Q;
    AICc(i) = logprob + 2*Q + 2*Q*(Q+1)/(N-Q-1);
    BIC(i) =  logprob + Q*log(N);

end

% Identify optimal models based on selection criteria
%-------------------------------------------------------------------------------
optima = zeros(nMethods,1);
functionals = cell(nMethods,1);
for i = 1:nMethods
    switch methods{i}
        case 'aic'
            functional = AIC;
        case 'aicc'
            functional = AICc;
        case 'bic'
            functional = BIC;
        case 'rmsd'
            functional = RMSD;
    end
    [~,optimum] = min(functional);
    functionals{i} = functional(:);
    optima(i) = optimum;
end

if nMethods==1
    functionals = functionals{1};
end

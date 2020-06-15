%
% FITPARAMODEL Fits a time- or distance-domain parametric model to one (or several) signals
%
%   [param,fit,paramuq,modeluq,stats] = FITPARAMODEL(V,@model,t,par0,lb,ub)
%   __ = FITPARAMODEL(V,@model,r,K,par0,lb,ub)
%   __ = FITPARAMODEL(V,@model,r,K,par0)
%   __ = FITPARAMODEL(V,@model,r,K)
%   __ = FITPARAMODEL(V,@model,t,par0,lb,ub)
%   __ = FITPARAMODEL(V,@model,t,par0)
%   __ = FITPARAMODEL(V,@model,t)
%   __ = FITPARAMODEL({V1,V2,___},@model,{t1,t2,___},par0,lb,ub)
%   __ = FITPARAMODEL({V1,V2,___},@model,r,{K1,K2,___},par0,lb,ub)
%   __ = FITPARAMODEL(___,'Property',Values)
%
%   Fits the N-point signal (V) to a M-point parametric model (@model) given an
%   M-point distance/time axis (r/t). For distance-domain fitting, provide
%   the NxM point kernel matrix (K). The fitted model corresponds to a parametric model
%   calculated by the passed function handle (@model). 

%   A structure containing different statistical estimators of goodness of 
%   fit is returned as (stats).
%
%   The initial guess of the model parameters can be passed as an
%   argument (param0). This is optional for DeerLab model functions. If (@model)
%   is a user-defined function handle, (param0) is required. Similarly, the
%   upper and lower bounds for the parameters can be passed as the last
%   arguments (ub) and (lb). This is optional for all model functions.
%
%   Pass multiple signals/kernels to enable global fitting of a single parametric
%   model to all data. The global fit weights are automatically computed according
%   to their contribution to ill-posedness.
%
%   Additional (optional) arguments can be passed as name-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'Solver' - Solver to be used to solve the minimization problems
%           'lsqnonlin' - Non-linear constrained least-squares (toolbox)
%             'fmincon' - Non-linear constrained minimization (toolbox)
%             'nlsqbnd' - Non-linear constrained least-squares (free)
%       'fminsearchcon' - Non-linear constrained minimization (free)
%          'fminsearch' - Unconstrained minimization
%   'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                     global fitting.
%   'Algorithm' - Algorithm to be used by the solvers (see fmincon or
%                 lsqnonlin documentation)
%   'TolFun' - Optimizer function tolerance
%   'MaxIter' - Maximum number of optimizer iterations
%   'MaxFunEvals' - Maximum number of optimizer function evaluations
%   'MultiStart' - Number of starting points for global optimization
%   'ConfidenceLevel' - Confidence level(s) for confidence intervals
%   'Rescale' - Enable/Disable optimization of the signal scale
%   'Verbose' - Display options for the solvers:
%                 'off' - no information displayed
%                 'final' - display solver exit message
%                 'iter-detailed' - display state of solver at each iteration
%               See MATLAB doc optimoptions for detailed explanation
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [parfit,modelfit,parci,modelci,stats] = fitparamodel(V,model,ax,varargin)

% Input parsing & validation
%-------------------------------------------------------------------------------
if nargin<3
    error('At least three inputs (V,model,t) are required.');
end

% Parse input schemes
%-------------------------------------------------------------------------------
% Prepare empty containers
K = [];
par0 = [];
lb = [];
ub = [];

% Check if kernel is passed
Kpassed = false;
if ~isempty(varargin)
    input = varargin{1};
    if ~ischar(input)
        if all(size(input)>1) || iscell(input)
            Kpassed = true;
            K = varargin{1};
            varargin(1) = [];
        end
    end
end
isDistanceDomain = Kpassed;

% Parse the varargin cell array
optionstart = numel(varargin);
for i=1:numel(varargin)
    if ischar(varargin{i})
        optionstart = i-1;
        break
    end
    switch i
        case 1
            par0 = varargin{1};
        case 2
            lb = varargin{2};
        case 3
            ub = varargin{3};
    end
end
varargin(1:optionstart) = [];

% Validate required inputs
%-------------------------------------------------------------------------------
if ~iscell(V)
    V = {V(:)};
end
nSignals = numel(V);
for i = 1:nSignals
    if ~iscolumn(V{i})
        V{i} = V{i}(:);
    end
    if ~isreal(V{i})
        V{i} = real(V{i});
    end
    if ~isreal(V{i})
        error('Input signal cannot be complex.')
    end
    validateattributes(V{i},{'numeric'},{'nonempty'},mfilename,'V')
end

% Validate input axis
if ~iscell(ax)
    ax = {ax(:)};
end
nAxes = numel(ax);
for i = 1:nAxes
    if ~iscolumn(ax{i})
        ax{i} = ax{i}(:);
    end
    % If using distance domain, only one axis must be passed
    if ~isDistanceDomain
        % Is using time-domain,control that same amount of axes as signals are passed
        if numel(V{i})~=numel(ax{i})
            error('V and t arguments must fulfill numel(t)==numel(S).')
        end
    end
end

if isDistanceDomain
    % Input #5 fitparamodel(V,model,r,K,'Property',Value)
    if nargin>4 && ischar(par0)
        varargin = [{par0},varargin];
        par0 = [];
    end
end

% Check that parametric model is passed as function handle
if ~isa(model,'function_handle')
    error('Model must be a function handle.')
end

% Get information about the parametric model
%-------------------------------------------------------------------------------
try
    % Check whether model is a DeerLab model function
    paraminfo = model();
    if nargin(model)==2
        model = @(ax,param,idx) model(ax,param);
    else
        model = @(ax,param,idx) model(ax,param,idx);
    end
    
catch
    % If not, require user to pass the inital values
    if ~exist('par0','var') || isempty(par0) || ischar(par0)
        error('For this model, please provide the required inital guess parameters.')
    end
    % Wrap the function handle into a DeerLab model function
    model = paramodel(model,par0,[],[],isDistanceDomain);
    paraminfo = model();
end

% Validate kernel
%-------------------------------------------------------------------------------
if isDistanceDomain
    if ~iscell(K)
        K = {K};
    end
    for i = 1:nSignals
        if size(K{i},1)~=numel(V{i})
            error('The number of rows in K must match the number of elements in V.')
        end
    end
    if numel(K)~=nSignals
        error('The number of kernels and signals must be equal.')
    end
end

if nargin<5 || isempty(par0)
    % If user does not give parameters, use the defaults of the model
    par0 =  [paraminfo.Start];
elseif nargin>4 && ischar(par0)
    varargin = [{par0} varargin];
    par0 = [paraminfo.Start];
else
    validateattributes(par0,{'numeric'},{'2d','nonempty'},mfilename,'StartParameters')
end

% Parse the optional parameters in varargin
%-------------------------------------------------------------------------------
[Solver,Algorithm,maxIter,Verbose,maxFunEvals,TolFun,GlobalWeights,...
    MultiStart,Rescale] = parseoptions(varargin);

% Validate optional inputs
if isempty(MultiStart)
    MultiStart = 1;
else
    validateattributes(MultiStart,{'numeric'},{'scalar','nonnegative'},mfilename,'MultiStarts')
end
if isempty(TolFun)
    TolFun = 1e-10;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonnegative'},mfilename,'TolFun')
end
if isempty(maxFunEvals)
    maxFunEvals = 5000;
else
    validateattributes(maxFunEvals,{'numeric'},{'scalar','nonnegative'},mfilename,'MaxFunEvals')
end
if isempty(maxIter)
    maxIter = 3000;
else
    validateattributes(maxIter,{'numeric'},{'scalar','nonnegative'},mfilename,'MaxIter')
end
if isempty(Verbose)
    Verbose = 'off';
else
    validateattributes(Verbose,{'char'},{'nonempty'},mfilename,'Verbose')
end
if isempty(Rescale)
    Rescale = true;
else
    validateattributes(Rescale,{'logical'},{'scalar'},mfilename,'Rescale')
end

% Solver
OptimizationToolboxInstalled = optimtoolbox_installed();
if OptimizationToolboxInstalled
    DefaultSolver = 'lsqnonlin';
    SolverList = {'lsqnonlin','nlsqbnd','lmlsqnonlin'};
else
    DefaultSolver = 'lmlsqnonlin';
    SolverList = {'nlsqbnd','lmlsqnonlin'};
end
if isempty(Solver)
    Solver = DefaultSolver;
end
validateattributes(Solver,{'char'},{'nonempty'},mfilename,'Solver')
validatestring(Solver,SolverList);
if ~ispc && strcmp(Solver,'nlsqbnd')
    error('The ''nlsqbnd'' solver is only available for Windows systems.')
end

% Algorithm
if isempty(Algorithm)
    if strcmp(Solver,'lsqnonlin')
        Algorithm = 'trust-region-reflective';
    else
        Algorithm = 'interior-point';
    end
else
    validInputs = {'levenberg-marquardt','interior-point','trust-region-reflective','active-set','sqp'};
    Algorithm = validatestring(Algorithm,validInputs);
end

% Global weights
if ~isempty(GlobalWeights)
    validateattributes(GlobalWeights,{'numeric'},{'nonnegative'})
    if numel(GlobalWeights)~=nSignals
        error('The number of global fit weights and signals must be equal.')
    end
    GlobalWeights = GlobalWeights/sum(GlobalWeights);
end


% Preparation of objective functions, parameter ranges, etc
%-------------------------------------------------------------------------------

% Parse errors in the model function, and reformat them
% model = @(ax,Parameters,idx)errorhandler(model,'modelfcn',ax,Parameters,idx);

% Get weights of different signals for global fitting
if isempty(GlobalWeights)
    Weights = globalweights(V);
else
    Weights = GlobalWeights;
end
Weights = Weights(:).';

Labels = num2cell(1:nSignals);

% Prepare upper/lower bounds on parameter search range
if isempty(lb)
    lb = [paraminfo.Lower];
end
if isempty(ub)
    ub = [paraminfo.Upper];
end
if any(ub==realmax) || any(lb==-realmax)
    unboundedparams = true;
    warning('Some model parameters are unbounded. Use ''Lower'' and ''Upper'' options to pass parameter boundaries.')
else
    unboundedparams = false;
end
if  numel(par0)~=numel(ub) || ...
        numel(par0)~=numel(lb)
    error('The inital guess and upper/lower boundaries must have equal length.')
end
if any(ub<lb)
    error('Lower bound values cannot be larger than upper bound values.')
end

% Preprare multiple start global optimization if requested
if MultiStart>1 && unboundedparams
    error('Multistart optimization cannot be used with unconstrained parameters.')
end
MultiStartParameters = multistarts(MultiStart,par0,lb,ub);



% Configure solver
%-------------------------------------------------------------------------------
switch Solver
    
    case 'lsqnonlin'
        solverOpts = optimoptions(@lsqnonlin,'Algorithm',Algorithm,...
            'Display',Verbose,'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
            'TolFun',TolFun,'DiffMinChange',0,'DiffMaxChange',Inf);
        solverFcn = @lsqnonlin;
        
    case 'lmlsqnonlin'
        
        solverOpts = struct('Display',Verbose,'MaxIter',maxIter,'MaxFunEvals',...
            maxFunEvals,'TolFun',TolFun);
        solverFcn = @lmlsqnonlin;
        
    case 'nlsqbnd'
        
        solverOpts = optimset('Algorithm',Algorithm,'Display',Verbose,...
            'MaxIter',maxIter,'MaxFunEvals',maxFunEvals,...
            'TolFun',TolFun,'TolCon',1e-20,...
            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        solverFcn = @nlsqbnd;
end

% Run least-squares fitting
%-------------------------------------------------------------------------------
% Disable ill-conditioned matrix warnings
warning('off','MATLAB:nearlySingularMatrix')
warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:rankDeficientMatrix')

fvals = zeros(1,MultiStart);
parfits = cell(1,MultiStart);
for runIdx = 1:MultiStart
    
    par0 = MultiStartParameters(runIdx,:);
    
    [parfit,fval,~,exitflag]  = solverFcn(@ResidualsFcn,par0,lb,ub,solverOpts);
    if exitflag==0
        % if maxIter exceeded, doube iterations and continue
        solverOpts.MaxIter = 2*maxIter;
        solverOpts.MaxFunEvals = 2*maxFunEvals;
        [parfit,fval]  = solverFcn(@ResidualsFcn,parfit,lb,ub,solverOpts);
    end
    
    fvals(runIdx) = fval;
    parfits{runIdx} = parfit;
    
end

% Find global minimum from multiple runs
[~,globmin] = min(fvals);
parfit = parfits{globmin};

% Issue warnings if fitted parameter values are at search range boundaries
%-------------------------------------------------------------------------------
tol = 1e-5;
atLower = abs(parfit-lb)<tol;
atUpper = abs(parfit-ub)<tol;
fixedPar = lb==ub;
for p = 1:numel(parfit)
    % Skip if parameter is fixed
    if fixedPar(p), continue; end
    
    if atLower(p)
        warning('The fitted value of parameter %d (%s) is at the lower bound of the range.',p,paraminfo(p).Parameter);
    end
    if atUpper(p)
        warning('The fitted value of parameter %d (%s) is at the upper bound of the range.',p,paraminfo(p).Parameter);
    end
end

% Calculate parameter confidence intervals
%-------------------------------------------------------------------------------
calcParamUncertainty = nargout>2;
if calcParamUncertainty
    
    % Compute residual vector and estimate variance from that
    residuals = ResidualsFcn(parfit);
    sigma2 = var(residuals);
    
    % Calculate numerical estimates of the Jacobian and Hessian
    % of the negative log-likelihood
    jacobian = jacobianest(@ResidualsFcn,parfit);
    hessian = jacobian.'*jacobian;
    
    % Estimate the covariance matrix by means of the inverse of Fisher information matrix
    lastwarn('');
    covmatrix = sigma2.*inv(hessian);
    % Detect if there was a 'nearly singular' warning
    [~, warnId] = lastwarn;
    if strcmp(warnId,'MATLAB:nearlySingularMatrix') || strcmp(warnId,'MATLAB:singularMatrix')
        covmatrix = sigma2.*sparse(pinv(full(hessian)));
        lastwarn('');
    end
    
    % Construct confidence interval structure
    parci = uqst('covariance',parfit,covmatrix,lb,ub);
    
end

% Evaluate fitted model and model confidence intervals
%-------------------------------------------------------------------------------
computeFittedModel = nargout>1;
if computeFittedModel
    modelfit = cell(nAxes,1);
    for i = 1:nAxes
        modelfit{i} = model(ax{i},parfit,Labels{i});
        if Rescale
            if isDistanceDomain
                scale = (K{i}*modelfit{i})\V{i};
            else
                scale = modelfit{i}\V{i};
            end
            modelfit{i} = scale*modelfit{i};
        end
    end
end

computeModelCI = nargout>3;
if computeModelCI
    modelci = cell(nAxes,1);
    %Loop over different signals
    for i = 1:nAxes
        if isDistanceDomain
            lb = zeros(numel(ax{i}),1);
        else
            lb = zeros(numel(ax{i}),1) - realmax;
        end
        ub = realmax + zeros(numel(ax{i}),1);
        modelci{i} = parci.propagate(@(par)model(ax{i},par,Labels{i}),lb,ub);
    end
end

% Calculate goodness of fit
%-------------------------------------------------------------------------------
computeStats = nargout>4;
if computeStats
    stats = cell(nSignals,1);
    for i = 1:nSignals
        if isDistanceDomain
            Vfit = K{i}*modelfit{1};
        else
            Vfit = modelfit{i};
        end
        Ndof = numel(V{i}) - numel(par0);
        stats{i} = gof(V{i},Vfit,Ndof);
    end
    
    if nSignals==1
        stats = stats{1};
    end
end

if nAxes==1
    if computeFittedModel
        modelfit = modelfit{1};
    end
    if computeModelCI
        modelci = modelci{1};
    end
end

% Set the warnings back on
warning('on','MATLAB:nearlySingularMatrix')
warning('on','MATLAB:singularMatrix')
warning('on','MATLAB:rankDeficientMatrix')

% Function that provides vector of residuals, which is the objective
% function for the least-squares solvers
    function r = ResidualsFcn(p)
        r_ = cell(nSignals,1);
        t = ax{1};
        for iSignal = 1:nSignals
            if nAxes>1
                t = ax{iSignal};
            end
            if isDistanceDomain
                sim = K{iSignal}*model(t,p,iSignal);
            else
                sim = model(t,p,iSignal);
            end
            if Rescale
                sim = (sim\V{iSignal})*sim;
            end
            r_{iSignal} = Weights(iSignal)*(V{iSignal}-sim);
        end
        r = vertcat(r_{:});
    end

end

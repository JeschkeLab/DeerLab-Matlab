%
%  SELREGPARAM Selection of optimal regularization parameter
%
%    alphaopt = SELREGPARAM(V,K,r,regtype,method)
%    alphaopt = SELREGPARAM(V,K,r,regtype,methodlist)
%    [alphaopt,F,alphas] = SELREGPARAM(___)
%    ___ = SELREGPARAM({V1,V2,...},{K1,K2,...},r,___)
%    alpha = SELREGPARAM(___,Name,Value)
%
%  Returns the optimal regularization parameter (alphaopt) from a range of
%  regularization parameter candidates. The parameter for the regularization
%  type given by (type) is computed based on the input signal (V) and the
%  dipolar kernel (K). The method employed for the selection of the 
%  regularization parameter is given in (method).
%
%  Inputs:
%    V          signal
%    K          dipolar kernel matrix
%    regtype    regularization type: 'tikhonov', 'tv', 'huber'
%    method     selection method: 'lr','lc','cv','gcv','rgcv','srgcv','aic',...
%                  'bic','aicc','rm','ee','ncp','gml','mcl'
%    methodlist cell array of selection methods, e.g. {'aic','bic','gcv'}
%   
%  Outputs:
%    alphaopt   optimal alpha, or alphas, determined
%    F          list of evaluated model-selection functionals
%    alphas     vector of evaluated alpha values
%
%    If multiple selection methods are passed as a cell array of strings,
%    the function returns (alpha) as an N-point array of optimal
%    regularization parameters corresponding to the input methods,
%
%    Passing multiple signals/kernels enables selection of the regularization
%    parameter for global fitting of the regularization model to a
%    single distribution. The global fit weights are automatically computed
%    according to their contribution to ill-posedness.
%
%  Name-value pairs:
%
%    'NonNegConstrained' - True/false to enforce non-negativity (default=true)
%    'RegOrder' - Order of the regularization operator L (default = 2).
%    'Refine' - True/false to enforce a second search around the optimal value
%               with a finer grid to achieve a better value of the optimum.
%    'HuberParameter' - Huber parameter used in the 'huber' model (default = 1.35).
%    'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                      global fitting regularization.
%    'TolFun' - Optimizer function tolerance.
%    'NoiseLevel' - Array of noise levels of the input signals.
%    'Range' - Range of alpha-value candidates to evaluate
%    'Search' - Search method to use: 'fminbnd' (default), 'golden', 'grid'
%

%  This file is a part of DeerLab. License is MIT (see LICENSE.md).
%  Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [alphaOpt,Functionals,alphaRanges,Residuals,Penalties] = selregparam(V,K,r,RegType,SelectionMethod,varargin)

%  Parse & validate required input
%-------------------------------------------------------------------------------

if nargin<4
    error('selregparam needs 5 inputs: V, K, RegType, SelectionMethod.')
end

% Turn off warnings to avoid ill-conditioned warnings
warning('off','MATLAB:nearlySingularMatrix')

% Check if user requested some options via name-value input
[TolFun,NonNegConstrained,NoiseLevel,GlobalWeights,HuberParameter,alphaRange,RegOrder,SearchMethod] ...
    = parseoptions(varargin);

if ~iscell(V)
    V = {V};
end
if ~iscell(K)
    K = {K};
end
if numel(K)~=numel(V)
    error('The number of kernels and signals must be equal.')
end

% Validate input signals
for i = 1:length(V)
    if ~iscolumn(V{i})
        V{i} = V{i}.';
    end
    if ~isreal(V{i})
        V{i} = error('Input signal(s) cannot be complex.');
    end
    if length(V{i})~=size(K{i},1)
        error('K and signal arguments must fulfill size(K,1)==length(S).')
    end
    validateattributes(V{i},{'numeric'},{'nonempty'},mfilename,'S')
    validateattributes(K{i},{'numeric'},{'nonempty'},mfilename,'K')
end

% Validate the selection methods input
allMethods = false;
allowedMethodInputs = {'lr','lc','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
if iscell(SelectionMethod)
    for i=1:length(SelectionMethod)
        if strcmp(SelectionMethod{i},'all')
            SelectionMethod = allowedMethodInputs;
            allMethods = true;
            break;
        end
        validateattributes(SelectionMethod{i},{'char'},{'nonempty'})
        SelectionMethod{i} = validatestring(SelectionMethod{i},allowedMethodInputs);
    end
else
    validateattributes(SelectionMethod,{'char'},{'nonempty'})
    if strcmp(SelectionMethod,'all')
        SelectionMethod = allowedMethodInputs;
        allMethods = true;
    else
        SelectionMethod = validatestring(SelectionMethod,allowedMethodInputs);
        SelectionMethod = {SelectionMethod};
    end
end
LcurveMethods = any(strcmpi(SelectionMethod,'lr')) || any(strcmpi(SelectionMethod,'lc'));

% Validate global weights
if ~isempty(GlobalWeights)
    validateattributes(GlobalWeights,{'numeric'},{'nonnegative'})
    if numel(GlobalWeights) ~= numel(V)
        error('The same number of global fit weights as signals must be passed.')
    end
    % Normalize weights
    GlobalWeights = GlobalWeights/sum(GlobalWeights);
else
    GlobalWeights = globalweights(V);
end

% Parse & validate optional input
%-------------------------------------------------------------------------------
warning('off','all')

if isempty(RegOrder)
    RegOrder = 2;
else
    validateattributes(RegOrder,{'numeric'},{'scalar','nonnegative'})
end

% Validate TolFun input
if isempty(TolFun)
    TolFun = 1e-9;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'nonNegLSQsolTol')
end

% Validate RegType input
if isempty(RegType)
    RegType = 'tikhonov';
else
    validateattributes(RegType,{'char'},{'nonempty'})
    allowedInput = {'tikhonov','tv','huber'};
    RegType = validatestring(RegType,allowedInput);
end
if isempty(HuberParameter)
    HuberParameter = 1.35;
else
    validateattributes(HuberParameter,{'numeric'},{'scalar','nonempty','nonnegative'})
end

% Validate NoiseLevel input
if isempty(NoiseLevel)
    for i = 1:length(V)
        NoiseLevel(i) = noiselevel(V{i});
    end
else
    if length(NoiseLevel)~=length(V)
        error('The same number of noise levels as signals must be provided.')
    end
    validateattributes(NoiseLevel,{'numeric'},{'scalar','nonempty','nonnegative'})
end

% Validate NonNegConstrained input
if isempty(NonNegConstrained)
    NonNegConstrained = true;
else
    validateattributes(NonNegConstrained,{'logical'},{'nonempty'},mfilename,'NonNegConstrained')
end

% Validate Search input
if isempty(SearchMethod)
    if LcurveMethods
        SearchMethod = 'grid';
    else
        SearchMethod = 'fminbnd';
    end
else
    validateattributes(SearchMethod,{'char'},{'nonempty'},mfilename,'Search')
    SearchMethod = validatestring(SearchMethod,{'grid','fminbnd'});
end

if allMethods && ~strcmp(SearchMethod,'grid')
    LcurveMethods = false;
    SelectionMethod(strcmp(SelectionMethod,'lr')) = [];
    SelectionMethod(strcmp(SelectionMethod,'lc')) = [];
end

if LcurveMethods && strcmp(SearchMethod,'fminbnd')
    error('The ''lr'' and ''lc'' selection methods work only with selregparam(...,''Search'',''grid'').')
end


%  Preparations
%-------------------------------------------------------------------------------
nr = size(K{1},2);
% Get regularization operator
L = regoperator(r,RegOrder);
% Get range of potential alpha values candidates
if isempty(alphaRange)
    alphaRange = regparamrange(K{1},L); %Does not work well with non-uniform r vectors
%     alphaRange = logspace(-10,10,100);
else
    validateattributes(alphaRange,{'numeric'},{'nonempty','nonnegative'},mfilename,'RegParamRange')
end

% Evaluate functional over search range, using specified search method
%-------------------------------------------------------------------------------
Functionals = cell(1,numel(SelectionMethod));
Residuals = cell(1,numel(SelectionMethod));
Penalties = cell(1,numel(SelectionMethod));
alphaRanges = cell(1,numel(SelectionMethod));
alphaOpt = zeros(1,numel(SelectionMethod));
switch lower(SearchMethod)
    
    case 'fminbnd'
        
        lga_min = log10(min(alphaRange));
        lga_max = log10(max(alphaRange));
        lga_opt = zeros(numel(SelectionMethod),1);
        for m = 1:numel(SelectionMethod)
            fun = @(lga) evalalpha(10.^lga,SelectionMethod(m));
            % Optimize alpha via the fminbnd function
            [lga_opt(m),~,history] = fminbnd_alpha(fun,lga_min,lga_max);
            % Get the vectors of evaluated alphas and functional values
            Functionals{m} = history(:,2);
            alphaRanges{m} = 10.^history(:,1);
            % Get the residual and penalty values 
            [~,Residuals{m},Penalties{m}] = fun(lga_opt(m));
        end
        alphaOpt = 10.^lga_opt;
    
    case 'grid'
        %-----------------------------------------------------------------------
        %  Grid search
        %-----------------------------------------------------------------------
        
        for a = 1:numel(alphaRange)
            [Functional(a,:),Residual(a,:),Penalty(a,:)] = evalalpha(alphaRange(a),SelectionMethod);
        end
        
        %  Grid-search specific selection methods
        for m = 1:length(SelectionMethod)
            for s = 1:length(V)
                switch lower(SelectionMethod{m})
                    
                    case 'lr' % L-curve minimum-radius method (LR)
                        Eta = log(Penalty(:,s));
                        Rho = log(Residual(:,s));
                        dd = @(x)(x-min(x))/(max(x)-min(x));
                        functional_ = dd(Rho).^2 + dd(Eta).^2;
                        
                    case 'lc' % L-curve maximum-curvature method (LC)
                        d1Residual = gradient(log(Residual(:,s)));
                        d2Residual = gradient(d1Residual);
                        d1Penalty = gradient(log(Penalty(:,s)));
                        d2Penalty = gradient(d1Penalty);
                        functional_ = (d1Residual.*d2Penalty - d2Residual.*d1Penalty)./(d1Residual.^2 + d1Penalty.^2).^(3/2);
                        
                    otherwise
                        functional_ = 0;
                end
                Functional(:,m) = Functional(:,m) + GlobalWeights(s)*functional_;
            end
        end
        
        for m = 1:numel(SelectionMethod)
            % Find index of selection functional minimum
            [~,idx] = min(Functional(:,m));
            % Store the corresponding regularization parameter
            alphaOpt(m) = alphaRange(idx);
            Functionals{m} = Functional(:,m);
            alphaRanges{m} = alphaRange(:);
            Residuals{m} = sum(Residual,2);
            Penalties{m} = sum(Penalty,2);
        end
end

% Do not return cell array if only one method is used
if numel(Functionals)==1
    Functionals = Functionals{1};
    Residuals = Residuals{1};
    alphaRanges = alphaRanges{1};
    Penalties = Penalties{1};
    Functionals = Functionals(:);
    Residuals = Residuals(:);
    Penalties = Penalties(:);
end


% Turn warnings back on
warning('on','MATLAB:nearlySingularMatrix');


%--------------------------------------------------------------------------
    function [Functional,Residual,Penalty] = evalalpha(alpha,SelectionMethod)
                
        %-----------------------------------------------------------------------
        %  Pseudo-Inverses and Ps
        %-----------------------------------------------------------------------
        
        [KtKreg,KtV] = lsqcomponents(V,K,L,alpha,RegType,HuberParameter,GlobalWeights);
        if NonNegConstrained
            P = fnnls(KtKreg,KtV,[],TolFun);
        else
            P = KtKreg\KtV;
        end
        for idx = 1:length(V)
            KtKreg = lsqcomponents(V{idx},K{idx},L,alpha,'tikhonov',HuberParameter,GlobalWeights);
            PseudoInverse{idx} = KtKreg\K{idx}.';
            switch lower(RegType)
                case 'tikhonov'
                    Penalty(idx) = norm(L*P);
                case 'tv'
                    Penalty(idx) = sum(sqrt((L*P).^2 + 1e-24));
                case 'huber'
                    Penalty(idx) = sum(sqrt((L*P/HuberParameter).^2 + 1 ) - 1);
            end
            Residual(idx) = norm(K{idx}*P - V{idx});
            InfluenceMatrix{idx} = K{idx}*PseudoInverse{idx};
        end
        
        %-----------------------------------------------------------------------
        %  Selection methods for optimal regularization parameter
        %-----------------------------------------------------------------------
        
        % If multiple selection methods are requested then process them sequentially
        Functional = zeros(length(SelectionMethod),1);
        for im = 1:length(SelectionMethod)
            for is = 1:length(V)
                nr = length(V{is});
                switch lower(SelectionMethod{im})
                    
                    case 'cv' % Cross validation (CV)
                        InfluenceDiagonal = diag(InfluenceMatrix{is});
                        f_ = sum(abs(V{is} - K{is}*(P)./(ones(nr,1) - InfluenceDiagonal)).^2);
                        
                    case 'gcv' % Generalized Cross Validation (GCV)
                        f_ = Residual(is)^2/((1 - trace(InfluenceMatrix{is})/nr)^2);
                        
                    case 'rgcv' % Robust Generalized Cross Validation (rGCV)
                        TuningParameter = 0.9;
                        f_ = Residual(is)^2/((1 - trace(InfluenceMatrix{is})/nr)^2)*(TuningParameter + (1 - TuningParameter)*trace(InfluenceMatrix{is}^2)/nr);
                        
                    case 'srgcv' % Strong Robust Generalized Cross Validation (srGCV)
                        TuningParameter = 0.8;
                        f_ = Residual(is)^2/((1 - trace(InfluenceMatrix{is})/nr)^2)*(TuningParameter + (1 - TuningParameter)*trace(PseudoInverse{is}'*PseudoInverse{is})/nr);
                        
                    case 'aic' % Akaike information criterion (AIC)
                        Criterion = 2;
                        f_ = nr*log(Residual(is)^2/nr) + Criterion*trace(InfluenceMatrix{is});
                        
                    case 'bic' % Bayesian information criterion (BIC)
                        Criterion = log(nr);
                        f_ = nr*log(Residual(is)^2/nr) + Criterion*trace(InfluenceMatrix{is});
                        
                    case 'aicc' % Corrected Akaike information criterion (AICC)
                        Criterion = 2*nr/(nr-trace(InfluenceMatrix{is})-1);
                        f_ = nr*log(Residual(is)^2/nr) + Criterion*trace(InfluenceMatrix{is});
                        
                    case 'rm' % Residual method (RM)
                        Scaling = K{is}.'*(eye(size(InfluenceMatrix{is})) - InfluenceMatrix{is});
                        f_ = Residual(is)^2/sqrt(trace(Scaling'*Scaling));
                        
                    case 'ee' % Extrapolated Error (EE)
                        f_ = Residual(is)^2/norm(K{is}.'*(K{is}*P - V{is}));
                        
                    case 'ncp' % Normalized Cumulative Periodogram (NCP)
                        ResidualPeriodogram = abs(fft(K{is}*P - V{is})).^2;
                        WhiteNoisePowerSpectrum = zeros(length(ResidualPeriodogram),1);
                        ResidualPowerSpectrum = zeros(length(ResidualPeriodogram),1);
                        for j=1:length(ResidualPeriodogram) - 1
                            ResidualPowerSpectrum(j)  = norm(ResidualPeriodogram(2:j+1),1)/norm(ResidualPeriodogram(2:end),1);
                            WhiteNoisePowerSpectrum(j) = j/(length(ResidualPeriodogram) - 1);
                        end
                        f_ = norm(ResidualPowerSpectrum - WhiteNoisePowerSpectrum);
                        
                    case 'gml' % Generalized Maximum Likelihood (GML)
                        Treshold = 1e-9;
                        EigenValues = eig(eye(size(InfluenceMatrix{is})) - InfluenceMatrix{is});
                        EigenValues(EigenValues < Treshold) = 0;
                        NonZeroEigenvalues = real(EigenValues(EigenValues~=0));
                        f_ = V{is}'*(V{is} - K{is}*P)/nthroot(prod(NonZeroEigenvalues),length(NonZeroEigenvalues));
                        
                    case 'mcl' % Mallows' C_L (MCL)
                        f_ = Residual^2 + 2*NoiseLevel(is)^2*trace(InfluenceMatrix{is}) - 2*nr*NoiseLevel(is)^2;
                    
                    otherwise
                        f_ = 0;
                end
                Functional(im) = Functional(im) + GlobalWeights(is)*f_;
            end
        end
    end

end


% Wrapper for the fminbnd function
function [x, fval, history] = fminbnd_alpha(fun,lga_min,lga_max)
    history = [];
    options = optimset('OutputFcn', @recordhistory,'TolX',0.01);
    [x, fval] = fminbnd(fun,lga_min,lga_max,options);
        
    % At each iteration record the current state of the optimization
    function stop = recordhistory(x,optimvalues,~)
        stop = false;
        history = [history; [x optimvalues.fval]];
    end
end

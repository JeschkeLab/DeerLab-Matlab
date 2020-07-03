%
% FITMULTIMODEL Multi-component distributions SNNLS fit function
%
%   [Pfit,parfit,Puq,paruq,Nopt,fcnal,Peval,stats] = fitmultimodel(___)
%   [___] = fitmultimodel(V,K,r,@model,maxModels,method,lb,ub)
%   [___] = fitmultimodel(V,K,r,@model,maxModels,method,lb)
%   [___] = fitmultimodel(V,K,r,@model,maxModels,method)
%   [___] = fitmultimodel(V,Kmodel,r,@model,maxModels,method,lb,ub)
%   [___] = fitmultimodel(V,Kmodel,r,@model,maxModels,method,lb)
%   [___] = fitmultimodel(V,Kmodel,r,@model,maxModels,method)
%   [___] = fitmultimodel({V1,V2,___},{K1,K2,___},r,@model,maxModels,method,lb,ub)
%   [___] = fitmultimodel({V1,V2,___},Kmodel,r,@model,maxModels,method,lb,ub)
%   [___] = fitmultimodel(___,'Name',Value)
%
%   Fits a multi-model parametric distance distribution model to the dipolar
%   signal (V) using (@model) as the basis function, the dipolar kernel (K)
%   and distance axis (r). The function compares these multi-model distributions
%   with up to a maximum number of basis functions given by (Nmax) and
%   determines the optimum one using the model selection criterion given in
%   (method) ('aic', 'bic', or 'aicc'). 
%
%   Ex/     Pfit = fitmultimodel(V,K,r,@dd_gauss,5,'aic'); 
%
%   Additional kernel parameters can be fitted by specifying a kernel function handle
%   (Kmodel) which accepts an array of parameters and returns a dipolar kernel matrix.
%   These parameters will be fitted along the other parameters. 
%
%   By passing multiple signals, a global distance distribution will be fitted 
%   to all signals. If a kernel model is specified, the kernel function must 
%   return a cell array of kernel matrices. Otherwise, a cell array of
%   kernel matrices {K1,K2,___} must be specified.
%
%   Inputs:
%       V           - dipolar signal, M-element array
%       K           - dipolar kernel, MxN-element matrix
%       Kmodel      - dipolar kernel model, function handle
%       r           - distance axis, N-element array
%       model       - model basis function, function handle
%       maxModels   - maximal number of components in model, scalar
%       method      - selection functional 
%                           'aic'  Akaike information criterion
%                           'aicc' corrected Akaike information criterion
%                           'bic'  Bayesian information criterion
%                           'rmsd' Root-mean squared deviation
%       lb          - Parameter lower bounds, 2-element cell array
%                           lb{1}  Kmodel-parameters lower bounds
%                           lb{2}  model-parameters lower bounds
%       ub          - Parameter upper bounds
%                           ub{1}  Kmodel-parameters upper bounds
%                           ub{2}  model-parameters upper bounds
%
%   Outputs: 
%       Pfit        - fitted distance distribution, N-element array
%       parfit      - fitted parameters, 3-element cell array
%                           parfit{1}  kernel parameters
%                           parfit{2}  distribution basis function parameters
%                           parfit{3}  amplitudes of the components
%       Puq         - distance distribution uncertainty quantification, struct
%       paruq       - parameters uncertainty quantification,struct
%       Nopt        - Optimal number of components, scalar
%       fcnal       - Values of the evaluated functonal specified by method
%       Peval       - Fitted distributions for all numbers of components
%       stats       - Goodness of fit, struct
%
%   Name-value pairs:
%      'GlobalWeights' - Array of weighting coefficients for the individual signals in
%                        global fitting (default = automatic).
%       'normP'        - true/false; whether to normalize Pfit (default = true)
%
%       See "help snlls" for a detailed list of other name-value pairs
%       accepted by the function.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [Pfit,fitparam,Puq,paruq,Nopt,fcnal,Peval,stats] = fitmultimodel(V,Kmodel,r,model,maxModels,method,varargin)

if nargin<6
    error('At least 6 inputs requried: fitmultimodel(V,Kmodel,r,model,maxModels,method)')
end

%Pre-allocate static workspace variables to share between subfunctions
[nKparam,lbK,ubK,nSignals] = deal([]);

% Parse inputs in the varargin
[lb,ub,options] = parseinputs(varargin);

% Validate all inputs
validateinputs()

% Default optional settings
GlobalWeights = globalweights(V);
normP = true;

% Parse and validate options passed by the user, if the user has specified
% any, this call will overwrite the defaults above
parsevalidate(options)

% Extract information about the model
info = model();
nparam = numel([info.Start]);
if isempty(lb0)
    lb0 = [info.Lower];
end
if isempty(ub0)
    ub0 = [info.Upper];
end
paramnames = {info.Parameter};

if any(cellfun(@(str)any(strcmp(str,{'Center','Location'})),paramnames))
    % If the center of the basis function is a parameter...
    idx = find((cellfun(@(str)any(strcmp(str,{'Center','Location'})),paramnames)));
    % ... limit it to the distance axis range (stabilizes parameter search)
    ub0(idx) = max(r);
    lb0(idx) = min(r);
end

% Pre-allocate containers
[Vfit,Pfit,plin_,pnonlin_,nlin_ub_,nlin_lb_,lin_ub_,lin_lb_] = deal(cell(maxModels,1));
[aic,aicc,bic,rmsd] = deal(zeros(maxModels,1));

% Loop over number of components in model
% =======================================
for Nmodels=1:maxModels
    
    % Prepare non-linear model with N-components
    % ===========================================
    Knonlin = @(par)nonlinmodel(par,Nmodels);
    
    % Box constraints for the model parameters (non-linear parameters)
    nlin_lb = repmat(lb0,1,Nmodels);
    nlin_ub = repmat(ub0,1,Nmodels);
    
    % Add the box constraints on the non-linear kernel parameters
    nlin_lb = [lbK nlin_lb];
    nlin_ub = [ubK nlin_ub];
    
    % Start values of non-linear parameters
    rng(1)
    par0 = (nlin_ub-nlin_lb).*rand(1,numel(nlin_lb)) + nlin_lb;
    
    % Box constraints for the components amplitudes (linear parameters)
    lin_lb = zeros(1,Nmodels); % Non-negativity constraint
    lin_ub = inf(1,Nmodels); % Unbounded
    
    % Separable non-linear least-squares (SNLLS) fit
    % ===============================================
    options_ = [options,'includePenalty',false]; % Do not use regularization
    [pnonlin,plin] = snlls(V,Knonlin,par0,nlin_lb,nlin_ub,lin_lb,lin_ub,options_{:});
    % Store the parameters for later
    pnonlin_{Nmodels} = pnonlin;
    plin_{Nmodels} = plin;
    nlin_ub_{Nmodels} = nlin_ub;
    nlin_lb_{Nmodels} = nlin_lb;
    lin_ub_{Nmodels} = lin_ub;
    lin_lb_{Nmodels} = lin_lb;
    
    % Get fitted kernel
    Kfit = nonlinmodel(pnonlin,Nmodels);
    
    % Get fitted signal
    Vfit{Nmodels} = cellfun(@(K) K*plin(:),Kfit,'UniformOutput',false);
    
    % Get fitted distribution
    Pfit{Nmodels} = Pmodel(pnonlin,plin);
    
    % Likelihood estimators
    % =====================
    [aic(Nmodels),aicc(Nmodels),bic(Nmodels),rmsd(Nmodels)] = logestimators(Vfit{Nmodels},plin,pnonlin);
    
end

Peval = Pfit;

% Select the optimal model
% ========================
fcnal = eval(method);
[~,Nopt] = min(fcnal);
Pfit = Pfit{Nopt};
Vfit = Vfit{Nopt};
pnonlin = pnonlin_{Nopt};
plin = plin_{Nopt};
nlin_lb = nlin_lb_{Nopt};
nlin_ub = nlin_ub_{Nopt};
lin_lb = lin_lb_{Nopt};
lin_ub = lin_ub_{Nopt};

% Package the fitted parameters
% =============================
fitparam{1} = pnonlin(1:nKparam); % Kernel parameters
fitparam{2} = pnonlin(nKparam+1:end); % Components parameters
fitparam{3} = plin; % Components amplitudes

% Uncertainty quantification analysis (if requested)
% ==================================================
if nargout>2
    [Puq,paruq] = uncertainty();
end

% Goodness of fit
% ===============
Ndof = nKparam + nparam + Nopt;
stats = cellfun(@(V,Vfit)gof(V,Vfit,Ndof),V(:),Vfit(:));

% If requested re-normalize the distribution
if normP
    Pnorm = trapz(r,Pfit);
    Pfit = Pfit/Pnorm;
    if nargout>2
        Puq.ci = @(p) Puq.ci(p)/Pnorm;
    end
end


% =========================================================================
% =========================================================================
% =========================================================================

% Non-linear augmented kernel model
% ------------------------------------------------------------------
    function Knonlin = nonlinmodel(par,Nmodels)
        
        K = Kmodel(par(1:nKparam));
        if ~iscell(K), K = {K}; end
        Knonlin = cell(nSignals,1);
        for iSignal = 1:nSignals
            subset = nKparam;
            for iModel = 1:Nmodels
                subset = subset(end)+1:subset(end)+nparam;
                % Get Gauss basis functions
                Pbasis = model(r,par(subset));
                
                % Combine all non-linear functions into one
                Knonlin{iSignal} = [Knonlin{iSignal} K{iSignal}*Pbasis];
            end
        end
    end

% Muli-component distribution model
% ------------------------------------------------------------------
% This function constructs the distance distribution from a set of
% non-linear and linear parameters given certain number of components and
% their basis function.
    function Pfit = Pmodel(nlinpar,linpar)
        subset = nKparam;
        Pfit = 0;
        for iModel = 1:numel(linpar)
            subset = subset(end)+1:subset(end)+nparam;
            % Get Gauss basis functions
            Pfit = Pfit + linpar(iModel)*model(r,nlinpar(subset));
        end
    end

% Log-Likelihood Estimators
% ------------------------------------------------------------------
% Computes the estimated likelihood of a multi-component model being the
% optimal choice.
    function [AIC,AICc,BIC,RMSD] = logestimators(Vfit,plin,pnonlin)
        [rmsd_,logprob] = deal(0);
        nParams = numel(pnonlin) + numel(plin);
        Q = nParams + 1;
        for iSignal = 1:nSignals
            N = numel(V{iSignal});
            SSR = sum((V{iSignal} - Vfit{iSignal}).^2);
            logprob = logprob + GlobalWeights(iSignal)*N*log(SSR/N);
            rmsd_ = rmsd_ + GlobalWeights(iSignal)*sqrt(1/N*SSR);
        end
        % Compute the estimators
        RMSD = rmsd_;
        AIC =  logprob + 2*Q;
        AICc = logprob + 2*Q + 2*Q*(Q+1)/(N-Q-1);
        BIC =  logprob + Q*log(N);
    end

% Uncertainty Quantificiation
% ------------------------------------------------------------------
% This function performs the covariance-based uncertaint analysis for the
% optimal multi-component model. It generates the uncertainty structure for
% the fit parameters and propagates it to the model.
    function [Puq,paramuq] = uncertainty()
        % Get full augmented residual vector
        %====================================
        Knonlin = @(par)nonlinmodel(par,Nopt);
        res = [];
        for j=1:nSignals
            res = [res; V{j} - cellselect(Knonlin(pnonlin),j)*plin(:)];
        end
        Jlin = [];
        Jnonlin = [];
        for j=1:nSignals
            Jnonlin = [Jnonlin; jacobianest(@(p)GlobalWeights(j)*cellselect(Knonlin(p),j)*plin(:),pnonlin)];
            Jlin = [Jlin; GlobalWeights(j)*cellselect(Knonlin(pnonlin),j)];
        end
        J = [Jnonlin Jlin];
        
        % Estimate the heteroscedasticity-consistent covariance matrix
        covmatrix = hccm(J,res,'HC1');
        
        % Construct uncertainty quantification structure for fitted parameters
        paramuq = uqst('covariance',[pnonlin(:); plin(:)],covmatrix,[nlin_lb(:); lin_lb(:)],[nlin_ub(:); lin_ub(:)]);
        
        P_subset = 1:nKparam+nparam*Nopt;
        amps_subset = P_subset(end)+1:P_subset(end)+Nopt;
        Puq = paramuq.propagate(@(p) Pmodel(p(P_subset),p(amps_subset)),zeros(numel(r),1));
    end

% Input parsing
% ------------------------------------------------------------------
% This function parses the different input schemes
    function  [lb,ub,options] = parseinputs(inputs)
        % Prepare empty containers
        [lb,ub] = deal([]);
        % Parse the varargin cell array
        optionstart = numel(inputs);
        for i=1:numel(inputs)
            if ischar(inputs{i})
                optionstart = i-1;
                break
            end
            switch i
                case 1
                    lb = inputs{i};
                case 2
                    ub = inputs{i};
            end
        end
        inputs(1:optionstart) = [];
        options = inputs;
    end

% Input validation
% ------------------------------------------------------------------
% This function checks that all inputs are valid
    function validateinputs()
        if ~iscell(V)
            V = {V};
        end
        nSignals = numel(V);
        for ii=1:nSignals
            V{ii} = V{ii}(:);
            validateattributes(V{ii},{'numeric'},{'vector','nonempty'});
        end
        % Basic validation of inputs
        validateattributes(maxModels,{'numeric'},{'scalar','integer','nonnegative','nonzero'});
        validateattributes(r,{'numeric'},{'vector','nonnegative'});
        validateattributes(model,{'function_handle'},{'nonempty'});
        validateattributes(method,{'char'},{'nonempty'});
        method = validatestring(method,{'aic','aicc','bic','rmsd'});
        r = r(:).';
        % Check kernel model
        if ~iscell(Kmodel)
            Kmodel = {Kmodel};
        end
        if isa(Kmodel{1},'function_handle')
            Kmodel = Kmodel{1};
            nKparam = 0;
            failing = true;
            while failing
                nKparam = nKparam + 1;
                try
                    Kmodel(rand(nKparam,1));
                    failing = false;
                catch
                    failing = true;
                end
            end
        else
            nKparam = 0;
            Kmodel = @(~)Kmodel;
        end
        
        % Parse boundaries
        if isempty(lb)
            lb = {[],[]};
        end
        if isempty(ub)
            ub = {[],[]};
        end
        lbK = lb{1};
        ubK = ub{1};
        lb0 = lb{2};
        ub0 = ub{2};
        if numel(lbK)~=nKparam || numel(ubK)~=nKparam
            error('The upper/lower bounds of the kernel parameters must be %i-element arrays',nKparam)
        end
    end

% Parsing and validation of options
% ------------------------------------------------------------------
% This function parses the name-value pairs or structures with options
% passes to the main function in the varargin. If the options pass, the
% default values are overwritten by the user-specified values.
    function parsevalidate(varargin)
        
        validoptions = {'normP','GlobalWeights'};
        
        % Parse options
        [normP_,GlobalWeights_] = parseoptions(validoptions,varargin);
        
        if ~isempty(normP_)
            validateattributes(normP_,{'logical'},{'scalar'})
            normP = normP_;
        end
        if ~isempty(GlobalWeights_)
            validateattributes(GlobalWeights_,{'numeric'},{'nonnegative','nonempty'})
            GlobalWeights = GlobalWeights_/sum(GlobalWeights_);
        end
    end

    function cell = cellselect(cell,idx)
        if ~iscell(cell)
            return
        end
        cell = cell{idx};
    end

end


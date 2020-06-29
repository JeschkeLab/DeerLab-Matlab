%
%  FITREGMODEL  Fits a non-parametric distance distribution to one (or several)
%               time-domain signals, using regularization
%
%   [Pfit,Pci,alpha,stats] = FITREGMODEL(V,K,r,regtype,method)
%   ___ = FITREGMODEL(V,K,r,regtype,method)
%   ___ = FITREGMODEL(V,K,r,regtype,alpha)
%   ___ = FITREGMODEL({V1,V2,___},{K1,K2,___},r,regtype,alpha)
%   ___ = FITREGMODEL(___,'Name',Values)
%
%   Regularization of the N-point signal (V) to a M-point distance
%   distribution (P) given a M-point distance axis (r) and NxM point kernel
%   (K). The regularization parameter (alpha) controls the regularization
%   properties. Instead of passing a numerial value for the regularization parameter
%   (alpha), the name of a selection method (method) can be passed and
%   the regularization parameter will be automatically selected by means
%   of the selregparam function.
%
%   The type of regularization employed in FITREGMODEL is set by the regtype
%   input argument. The regularization models implemented in FITREGMODEL are:
%          'tikhonov' -   Tikhonov regularization
%          'tv'       -   Total variation regularization
%          'huber'    -   pseudo-Huber regularization
%
%   Passing multiple signals/kernels enables global fitting of the
%   a single-distribution model to all data. The global fit weights
%   are automatically computed according to their contribution to ill-posedness.
%
%   Additional (optional) arguments can be passed as property-value pairs.
%
%  Input:
%    V        time-domain signal to fit (N-element vector)
%    K        dipolar kernel(NxM-element matrix)
%    r        distance axis, in nanometers (M-element vector)
%    regtype  regularization penalty type (string)
%                'tikhonov' -   Tikhonov regularization
%                'tv'       -   Total variation regularization
%                'huber'    -   pseudo-Huber regularization
%    method   regularization parameter selection method (string) (see help selregparam)
%    alpha    regularization parameter value (scalar)
%
%  Output:
%    Pfit     fitted distance distribution (M-element vector)
%    Pci      fitted distribution uncertainty quantification (struct)
%    alpha    regularization parameter used for the analysis
%    stats    goodness of fit statistical estimators
%
%  Name-value pairs:
%   'Solver'            - Solver to be used to solve the minimization problems
%                            'fnnls' - Fast non-negative least-squares
%                            'lsqnonneg' - Non-negative least-squares
%                            'fmincon' - Non-linear constrained minimization
%                            'bppnnls' -  Block principal pivoting non-negative least-squares solver
%   'NonNegConstrained' - Enable/disable non-negativity constraint (true/false)
%   'OBIR'              - Use Osher's Bregman iterated regularization algorithm for fitting
%   'HuberParam'        - Huber parameter used in the 'huber' model (default = 1.35)
%   'NormP'             - true/false; whether to normalize Pfit (default = true)
%   'GlobalWeights'     - Array of weighting coefficients for the individual signals in
%                         global fitting.
%   'RegOrder'          - Order of the regularization operator L (default = 2)
%   'ConfidenceLevel'   - Confidence lvel(s) for the confidence intervals (0-1, default=0.95)
%   'TolFun'            - Optimizer function tolerance
%   'MaxIter'           - Maximum number of optimizer iterations
%   'MaxFunEvals'       - Maximum number of optimizer function evaluations
%   'Verbose'           - Display options for the solvers:
%                           'off' - no information displayed
%                           'final' - display solver exit message
%                           'iter-detailed' - display state of solver at each iteration
%                         See MATLAB doc optimoptions for detailed explanation
%   See "help selregparam" for further options.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [P,Pci,alpha,stats] = fitregmodel(V,K,r,RegType,alpha,varargin)

% Required Input
% ===============

parsevalidate_inputs(nargin);

% Optional Input
% ===============

% Default optional settings
Verbose = 'off';
RegOrder = 2;
TolFun = 1e-9;
runOBIR = false;
MaxIter = 2e7;
HuberParam = 1.35;
NormP = true;
MaxFunEvals = 2e7;
NonNegConstrained = true;
GlobalWeights = globalweights(V);
if strcmp(RegType,'custom')
    Solver = 'fmincon';
else
    Solver = 'fnnls';
end

% Parse and validate options passed by the user, if the user has specified
% any, this call will overwrite the defaults above
parsevalidate_options(varargin)

% Decide to do uncertainty quantification
getConfidenceIntervals = nargout>1 & ~strcmp(RegType,'custom');


% Regularization
% ==================

% Turn off warnings to avoid ill-conditioned warnings
warning('off','MATLAB:nearlySingularMatrix')

% Comppute regularization matrix
L = regoperator(r,RegOrder);

nr = size(K{1},2);
P0 = zeros(nr,1);

% If unconstrained regularization is requested then solve analytically
if ~NonNegConstrained && ~strcmp(Solver,'fmincon')
    Solver = 'analytical';
end

% If using LSQ-based solvers then precompute the KtK and KtS input arguments
if ~strcmp(Solver,'fmincon') || getConfidenceIntervals
    [KtKreg,KtV] = lsqcomponents(V,K,L,alpha,RegType,HuberParam,GlobalWeights);
end

if runOBIR

    % Osher's Bregman iterated regularization
    % ==========================================    
    % Fit using Osher's Bregman iterated regularization algorithm
    P = obir(V,K,L,RegType,alpha,GlobalWeights,HuberParam,Solver);
    
else
    
    % Regularization
    % ==========================================
    % Solve the regularization functional minimization problem
    switch lower(Solver)
        
        case 'analytical'
            P = zeros(nr,1);
            for i = 1:length(V)
                PseudoInverse = KtKreg\K{i}.';
                P = P + GlobalWeights(i)*PseudoInverse*V{i};
            end
            
        case 'lsqnonneg'
            solverOpts = optimset('Display','off','TolX',TolFun);
            P = lsqnonneg(KtKreg,KtV,solverOpts);
            
        case 'fnnls'
            [P,~,~,flag] = fnnls(KtKreg,KtV,P0,TolFun,Verbose);
            % In some cases, fnnls may return negatives if tolerance is too high
            if flag==-1
                %... in those cases continue from current solution
                [P,~,~,flag] = fnnls(KtKreg,KtV,P,1e-20);
            end
            if flag==-2
                warning('FNNLS cannot solve the problem. Regularization parameter may be too large.')
            end
            
        case 'bppnnls'
            P = nnls_bpp(KtKreg,KtV,KtKreg\KtV);
            
        case 'fmincon'
            if NonNegConstrained
                NonNegConst = zeros(nr,1);
            else
                NonNegConst = [];
            end
            if ~strcmp(RegType,'custom')
                RegFunctional = regfunctional(RegType,V,L,K,alpha,HuberParam);
            else
                % Parse errors in the analyzed function, and reformat them
                RegFunctional = @(P)errorhandler(RegFunctional,'regfcn',P);
            end
            
            GradObj = ~strcmp(RegType,'custom');
            fminconOptions = optimoptions(@fmincon,'SpecifyObjectiveGradient',GradObj,'MaxFunEvals',MaxFunEvals,'Display',Verbose,'MaxIter',MaxIter);
            [P,~,exitflag] =  fmincon(RegFunctional,P0,[],[],[],[],NonNegConst,[],[],fminconOptions);
            % Check how optimization exited...
            if exitflag == 0
                %... if maxIter exceeded (flag=0) then double iterations and continue from where it stopped
                fminconOptions = optimoptions(fminconOptions,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
                P  = fmincon(RegFunctional,P,[],[],[],[],NonNegConst,[],[],fminconOptions);
            end
    end
    
end

% Uncertainty quantification
% ================================
if getConfidenceIntervals
    warning('off','MATLAB:nearlySingularMatrix')
    
    % Estimate the contribution to P variance from the different signals
    res = [];
    J = [];
    for ii=1:numel(V)
        % Construct the full augmented residual
        res = [res; GlobalWeights(ii)*(V{i} - K{i}*P)];
        J = [J; GlobalWeights(ii)*K{i}];
    end
    Jreg = alpha*L;
    % Augment residual and Jacobian with regularization term
    res = [res; Jreg*P];
    J = [J; Jreg];
    
    covmat = hccm(J,res,'HC1');
    % Construct confidence interval structure for P
    NonNegConst = zeros(nr,1);
    Pci = uqst('covariance',P,covmat,NonNegConst,[]);
    
end


% Re-normalization
% ================================
if NormP
    Pnorm = trapz(r,P);
    P = P/Pnorm;
    if getConfidenceIntervals
        Pci.ci = @(p)Pci.ci(p)/Pnorm;
    end
end

% Goodness of fit
% ================================
% If requested compute Goodness of Fit for all signals
if nargout>3
    stats = cell(numel(V),1);
    for i=1:numel(V)
        Vfit = K{i}*P;
        % Estimate degrees of freedom from the inlfuence matrix
        KtKreg = lsqcomponents(V{i},K{i},L,alpha,RegType,HuberParam);
        H = K{i}*(KtKreg\K{i}.');
        Ndof = numel(V{i}) - trace(H);
        stats{i} = gof(V{i},Vfit,Ndof);
    end
    if numel(V)==1
        stats = stats{1};
    end
end

% Turn warnings back on
warning('on','MATLAB:nearlySingularMatrix')


% Parsing and validation of requried inputs
% ------------------------------------------------------------------
% This function parses the inputs of the main function.
    function parsevalidate_inputs(nInputs)
        
        if nInputs<3
            error('Not enough input arguments.')
        end
        if nInputs<4 || isempty(RegType)
            RegType = 'tikhonov';
        elseif isa(RegType,'function_handle')
            RegFunctional = RegType;
            RegType = 'custom';
        else
            validateattributes(RegType,{'char'},{'nonempty'})
            allowedInput = {'tikhonov','tv','huber','custom'};
            RegType = validatestring(RegType,allowedInput);
        end
        if  nInputs<5 || isempty(alpha)
            alpha = 'aic';
        end
        if ~iscell(V)
            V = {V};
        end
        if ~iscell(K)
            K = {K};
        end
        if isa(alpha,'char')
            alpha = selregparam(V,K,r,RegType,alpha,varargin);
        else
            validateattributes(alpha,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
        end
        validateattributes(r,{'numeric'},{'nonempty','nonnegative'},mfilename,'r')
        
        if numel(K)~=numel(V)
            error('The number of kernels must be equal to the number of kernels.')
        end
        for j = 1:numel(V)
            if ~iscolumn(V{j})
                V{j} = V{j}.';
            end
            if ~isreal(V{j})
                V{j} = real(V{j});
            end
            if length(V{j})~=size(K{j},1)
                error('The number of rows in K must match the number of elements in V.')
            end
            validateattributes(V{j},{'numeric'},{'nonempty'},mfilename,'S')
        end
    end

% Parsing and validation of options
% ------------------------------------------------------------------
% This function parses the name-value pairs or structures with options
% passes to the main function in the varargin. If the options pass, the
% default values are overwritten by the user-specified values.
    function parsevalidate_options(options)
        
        optionalProperties = {'TolFun','Solver','NonNegConstrained','Verbose',...
            'MaxFunEvals','MaxIter','HuberParam','GlobalWeights',...
            'OBIR','RegOrder','NormP'};
        [TolFun_,Solver_,NonNegConstrained_,Verbose_,MaxFunEvals_,MaxIter_,HuberParam_,...
            GlobalWeights_,runOBIR_,RegOrder_,NormP_] = parseoptions(optionalProperties,options);
        
        if ~isempty(GlobalWeights_)
            validateattributes(GlobalWeights_,{'numeric'},{'nonnegative'})
            if numel(GlobalWeights_) ~= numel(V)
                error('The same number of global fit weights as signals must be passed.')
            end
            % Normalize weights
            GlobalWeights = GlobalWeights_/sum(GlobalWeights_);
        end
        
        if ~isempty(Verbose_)
            validateattributes(Verbose_,{'char'},{'nonempty'},mfilename,'Verbose')
            Verbose = Verbose_;
        end
        if ~isempty(RegOrder_)
            validateattributes(RegOrder_,{'numeric'},{'scalar','nonnegative'})
            RegOrder = RegOrder_;
        end
        if ~isempty(TolFun_)
            validateattributes(TolFun_,{'numeric'},{'scalar','nonempty','nonnegative'},'fitregmodel','TolFun')
            TolFun = TolFun_;
        end
        if ~isempty(Solver_)
            validateattributes(Solver_,{'char'},{'nonempty'})
            allowedInput_ = {'analytical','fnnls','lsqnonneg','bppnnls','fmincon'};
            Solver = validatestring(Solver_,allowedInput_);
        end
        
        if ~isempty(runOBIR_)
            validateattributes(runOBIR_,{'logical'},{'nonempty'},'fitregmodel','runOBIR')
            runOBIR = runOBIR_;
        end
        
        if ~isempty(MaxIter_)
            validateattributes(MaxIter_,{'numeric'},{'scalar','nonempty'},mfilename,'MaxIter')
            MaxIter = MaxIter_;
        end
        
        if ~isempty(HuberParam_)
            validateattributes(HuberParam_,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'MaxFunEvals')
            HuberParam = HuberParam_;
        end
        
        if ~isempty(NormP_)
            validateattributes(NormP_,{'logical'},{'scalar'},mfilename,'NormP');
            NormP = NormP_;
        end
        
        if ~isempty(MaxFunEvals_)
            validateattributes(MaxFunEvals_,{'numeric'},{'scalar','nonempty'},mfilename,'MaxFunEvals')
            MaxFunEvals = MaxFunEvals_;
        end
        
        if ~isempty(NonNegConstrained_)
            validateattributes(NonNegConstrained_,{'logical'},{'nonempty'},'fitregmodel','NonNegConstrained')
            NonNegConstrained = NonNegConstrained_;
        end
        
    end

end

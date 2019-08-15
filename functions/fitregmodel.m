% 
% FITREGMODEL Fits a distance distribution to one (or several) signals
%            by optimization of a regularization functional model.
%
%   P = FITREGMODEL(S,R,K,L,'type',ALPHA)
%   Regularization of the N-point signal S to a M-point distance
%   distribution P given a M-point distance axis R and NxM point kernel
%   K. The (M-2)xM point regularization matrix L and regularization
%   parameter ALPHA control the regularization properties.
%
%   The type of regularization employed in FITREGMODEL is set by the 'type'
%   input argument. The regularization models implemented in FITREGMODEL are:
%          'tikhonov' -   Tikhonov regularization
%          'tv'       -   Total variation regularization
%          'huber'    -   pseudo-Huber regularization
%
%   K = FITREGMODEL(...,'Property',Values)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order. 
%
%   'Solver' - Solver to be used to solve the minimization problems
%                      'fnnls' - Fast non-negative least-squares
%                      'lsqnonneg' - Non-negative least-squares 
%                      'fmincon' - Non-linear constrained minimization
%                      'bnnlps'
%
%   'NonNegConstrained' - Enable/disable non-negativity constraint (true/false)
%
%   'HuberParam' - Huber parameter used in the 'huber' model (default = 1.35).
%
%   'TolFun' - Optimizer function tolerance
%
%   'MaxIter' - Maximum number of optimizer iterations
%
%   'MaxFunEvals' - Maximum number of optimizer function evaluations   
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function Distribution = fitregmodel(Signal,DistanceAxis,Kernel,RegMatrix,RegType,RegParam,varargin)

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if nargin<5
    error('Not enough input arguments.')
end
if nargin<4 || isempty(RegType)
    RegType = 'tikhonov';
elseif isa(RegType,'function_handle')
    RegFunctional = RegType;
    RegType = 'custom';
else
    validateattributes(RegType,{'char'},{'nonempty'})
    allowedInput = {'tikhonov','tv','huber'};
    validatestring(RegType,allowedInput);
end

if strcmp(RegType,'custom')
    GradObj = false;
else
    GradObj = true;
end
validateattributes(RegMatrix,{'numeric'},{'nonempty','2d'},mfilename,'RegMatrix')
validateattributes(RegParam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
validateattributes(DistanceAxis,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'DistanceAxis')
if numel(unique(round(diff(DistanceAxis),12)))~=1
    error('Distance axis must be a monotonically increasing vector.')
end

%--------------------------------------------------------------------------
% Parse & Validate Optional Input
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[TolFun,Solver,NonNegConstrained,MaxFunEvals,MaxIter,HuberParam] = parseoptional({'TolFun','Solver','NonNegConstrained','MaxFunEvals','MaxIter','HuberParam'},varargin);

if isempty(TolFun)
    TolFun = 1e-9;
else
    validateattributes(TolFun,{'numeric'},{'scalar','nonempty','nonnegative'},'regularize','nonNegLSQsolTol')
end
if isempty(Solver)
    Solver = 'fnnls';
else
    validateattributes(Solver,{'char'},{'nonempty'})
    allowedInput = {'analytical','fnnls','lsqnonneg','bppnnls','fmincon'};
    validatestring(Solver,allowedInput);
end

if isempty(MaxIter)
    MaxIter = 20000000;
else
    validateattributes(MaxIter,{'numeric'},{'scalar','nonempty'},mfilename,'MaxIter')
end

if isempty(HuberParam)
    HuberParam = 1.35;
else
    validateattributes(HuberParam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'MaxFunEvals')
end

if isempty(MaxFunEvals)
    MaxFunEvals = 2000000;
else
    validateattributes(MaxFunEvals,{'numeric'},{'scalar','nonempty'},mfilename,'MaxFunEvals')
end

if isempty(NonNegConstrained)
    NonNegConstrained = true;
else
    validateattributes(NonNegConstrained,{'logical'},{'nonempty'},'regularize','NonNegConstrained')
end
if ~iscell(Signal)
    Signal = {Signal};
end
if ~iscell(Kernel)
    Kernel = {Kernel};
end
if length(Kernel)~=length(Signal)
    error('The number of kernels and signals must be equal.')
end
for i=1:length(Signal)
    if ~iscolumn(Signal{i})
        Signal{i} = Signal{i}.';
    end
    if ~isreal(Signal{i})
        Signal{i} = real(Signal{i});
    end
    if length(Signal{i})~=size(Kernel{i},1)
        error('Kernel and signal arguments must fulfill size(Kernel,1)==length(Signal).')
    end
    validateattributes(Signal{i},{'numeric'},{'nonempty'},mfilename,'Signal')
end

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end

hashKey = datahash({Signal,Kernel,RegMatrix,RegType,RegParam,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [Distribution] = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
%Regularization processing
%--------------------------------------------------------------------------

Dimension = length(RegMatrix);
InitialGuess = zeros(Dimension,1);
dr = mean(diff(DistanceAxis));

%If unconstrained regularization is requested then solve analytically
if ~NonNegConstrained && ~strcmp(Solver,'fmincon')
    Solver = 'analytical';
end

%If using LSQ-based solvers then precompute the KtK and KtS input arguments
if ~strcmp(Solver,'fmincon')
    [Q,KtS,weights] =  lsqcomponents(Signal,Kernel,RegMatrix,RegParam,RegType,HuberParam);
end

switch lower(Solver)
    
    case 'analytical'
        Distribution = zeros(Dimension,1);
        for i=1:length(Signal)
        PseudoInverse = Q\Kernel{i}.';
        Distribution = Distribution + weights(i)*PseudoInverse*Signal{i};
        end
    case 'lsqnonneg'
        solverOpts = optimset('Display','off','TolX',TolFun);
        Distribution = lsqnonneg(Q,KtS,solverOpts);

    case 'fnnls'
        Distribution = fnnls(Q,KtS,InitialGuess,TolFun);
        %In some cases, fnnls may return negatives if tolerance is to high
        if any(Distribution < 0)
            %... in those cases continue from current solution
            Distribution = fnnls(Q,KtS,Distribution,1e-20);
        end
    case 'bppnnls'
        Distribution = nnls_bpp(Q,KtS,Q\KtS);
        
    case 'fmincon'
        %Constrained Tikhonov/Total variation/Huber regularization
        if NonNegConstrained
            NonNegConst = zeros(Dimension,1);
        else
            NonNegConst = [];
        end
        if ~strcmp(RegType,'custom')
            RegFunctional = regfunctional(RegType,Signal,RegMatrix,Kernel,RegParam,HuberParam);
        end
        constraint = @(x)unityconstraint(x,dr);
        fminconOptions = optimoptions(@fmincon,'SpecifyObjectiveGradient',GradObj,'MaxFunEvals',MaxFunEvals,'Display','off','MaxIter',MaxIter);
        [Distribution,~,exitflag] =  fmincon(RegFunctional,InitialGuess,[],[],[],[],NonNegConst,[],constraint,fminconOptions);
        %Check how optimization exited...
        if exitflag == 0
            %... if maxIter exceeded (flag =0) then doube iterations and continue from where it stopped
            fminconOptions = optimoptions(fminconOptions,'MaxIter',2*MaxIter,'MaxFunEvals',2*MaxFunEvals);
            Distribution  = fmincon(ModelCost,Distribution,[],[],[],[],NonNegConst,[],constraint,fminconOptions);
        end
end

%Normalize distribution integral
Distribution = Distribution/sum(Distribution)/dr;

%Store output result in the cache
cachedData = addcache(cachedData,hashKey,Distribution);

%--------------------------------------------------------------------------
end
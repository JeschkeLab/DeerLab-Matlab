%
% OBIR Osher's Bregman-iterated regularization method
%
%   P = OBIR(V,K,r,'type',alpha)
%
%   OBIR of the N-point signal (V) to a M-point distance
%   distribution (P) given a M-point distance axis (r) and NxM point kernel
%   (K). The regularization parameter (alpha) controls the regularization
%   properties.
%
%   The type of regularization employed in OBIR is set by the 'type'
%   input argument. The regularization models implemented in OBIR are:
%       'tikhonov' -   Tikhonov regularization
%       'tv'       -   Total variation regularization
%       'huber'    -   pseudo-Huber regularization
%
%   P = OBIR(...,'Property',Value)
%   Additional (optional) arguments can be passed as name-value pairs.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function Pfit = obir(V,K,L,RegType,alpha,GlobalWeights,HuberParam,Solver)

NoiseLevelAim = zeros(1,numel(V));
for i=1:numel(V)
    NoiseLevelAim(i) = noiselevel(V{i});
end

Verbose = 'off';
TolFun = 1e-10;
MaxOuterIter = 5000;
DivergenceStop = false;

% Preparation
%-------------------------------------------------------------------------------

% Initialize
nr = size(L,2);
subGrad = zeros(nr,1);
Counter = 1;
Iteration = 1;
Pfit = zeros(nr,1);

% Osher's Bregman Iterated Algorithm
%-------------------------------------------------------------------------------
InitialGuess = zeros(nr,1);
diverged = zeros(numel(V),1);
semiconverged = zeros(numel(V),1);

% Precompute the KtK and KtS input arguments
[KtKreg,KtV] = lsqcomponents(V,K,L,alpha,RegType,HuberParam,GlobalWeights);

while Iteration <= MaxOuterIter
    
    % Store previous iteration distribution
    Pprev = Pfit;
    
    %Update
    KtV_ = KtV - subGrad;
    
    %Get solution of current Bregman iteration
    switch lower(Solver)
        case 'fnnls'
            Pfit = fnnls(KtKreg,KtV_,InitialGuess,TolFun,Verbose);
        case 'lsqnonneg'
            Pfit = lsqnonneg(KtKreg,KtV_,solverOpts);
        case 'bppnnls'
            Pfit = nnls_bpp(KtKreg,KtV_,KtKreg\KtV_);
    end
    
    % Update subgradient at current solution
    for i=1:numel(V)
        subGrad = subGrad + GlobalWeights(i)*K{i}.'*(K{i}*Pfit - V{i});
    end
    
    % Iteration control
    %--------------------------------------------------------------------------
    if Iteration == 1
        % If at first iteration, the residual deviation is already below the
        % noise deviation then impose oversmoothing and remain at first iteration
        if NoiseLevelAim  > std(K{1}*Pfit - V{1})
            alpha = alpha*2^Counter;
            Counter = Counter + 1;
            
            % Recompute the KtK and KtS input arguments with new alpha
            [KtKreg,KtV] = lsqcomponents(V,K,L,alpha,RegType,HuberParam);
        else
            % Once the residual deviation is above the threshold, then proceed
            % further with the Bregman iterations
            Iteration  = Iteration + 1;
        end
    else
        % For the rest of the Bregman iterations control the condition and stop
        % when fulfilled
        for i=1:numel(V)
         diverged(i) = std(K{i}*Pprev - V{i}) < std(K{i}*Pfit - V{i});
         semiconverged(i) = NoiseLevelAim(i) > std(K{i}*Pfit - V{i});
        end
        
        if all(semiconverged)
            break;
        else
            Iteration  = Iteration + 1;
        end
        % If residual deviation starts to diverge, stop
        if DivergenceStop && any(diverged)
            Pfit = Pprev;
            break;
        end
    end
    
end

end
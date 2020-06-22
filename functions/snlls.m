%
% SNNLS  Separable Non-linear Least Squares Solver
%
%   [pnlin,plin,paramuq] = SNNLS(y,Amodel,par0,lb,ub,lbl,ubl)
%   __ = SNNLS(y,Amodel,par0,lb,ub,lbl,ubl)
%   __ = SNNLS(y,Amodel,par0,lb,ub)
%   __ = SNNLS(y,Amodel,par0)
%   __ = SNNLS(y,Amodel,par0)
%   __ = SNNLS(___,'Property',Values,___)
%
%   Fits a linear set of parameters (x) and non-linear parameters (p) 
%   by solving the following non-linear least squares problem:
%           [x,p] = argmin || A(p)*x - y||^2
%                    s.t.   x in [lbl,ubl] 
%                           p in [lb,ub] 
%
%  Input:
%    y        N-element vector of input data to be fitted
%    Amodel   Function handle accepting non-linear parameters and returning
%             a NxM-element matrix 
%    par0     W-element vector, start values of the non-linear parameters
%    lb       W-element vector of lower bounds for the non-linear parameters
%    ub       W-element vector of upper bounds for the non-linear parameters
%    lbl      M-element vector of lower bounds for the linear parameters
%    ubl      M-element vector of lower bounds for the linear parameters
%
%  Output:
%    pnlin    fitted non-linear parameters
%    nlin     fitted linear parameters
%    paramuq  uncertainty quantification structure
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [nonlinfit,linfit,paramuq] = snlls(y,Amodel,par0,lb,ub,lbl,ubl,RegOrder)

if isempty(RegOrder)
    RegOrder = 2;
else
    validateattributes(RegOrder,{'numeric'},{'scalar','integer'},mfilename,'RegOrder')
end

if optimtoolbox_installed
    
    % Options for non-linear solver
    nonLinSolverFcn = @lsqnonlin;
    nonLinSolverOpts = optimoptions(@lsqnonlin,'Display','off','MaxIter',1e4,...
        'MaxFunEvals',1e4,'TolFun',1e-5,...
        'DiffMinChange',0,'DiffMaxChange',Inf);
    
    % Options for linear solver
    linSolverFcn = @(A,y,lbl,ubl,opts)lsqlin(A,y,[],[],[],[],lbl,ubl,[],opts);
    linSolverOpts = optimoptions(@lsqlin,'Display','off','MaxIter',1e4,'TolFun',1e-5);
else
    % Options for non-linear solver
    nonLinSolverFcn = @lmlsqnonlin;
    nonLinSolverOpts = struct('Display','off','MaxIter',1e4,...
        'MaxFunEvals',1e4,'TolFun',1e-5);
    
    % Options for linear solver
    linSolverFcn = @lsqlin_QP;
    linSolverOpts = [];
end

% Check if the linear problem is constrained
if ~isempty(lbl) || ~isempty(ubl)
    linearConstrained = true;
end

% Check for non-negativity constraints on the linear solution
if all(lbl==0) && isempty(ubl)
    nonNegativeOnly = true;
else
    nonNegativeOnly = false;
end

%Pre-allocate static workspace
illConditioned = [];
L = [];
regparam = 'aic';
par_prev = [];
regparam_prev = [];
alphaOptThreshold = 1e-3;
linfit = [];

% Run the non-linear solver
[nonlinfit,fval,~,exitflag]  = nonLinSolverFcn(@ResidualsFcn,par0,lb,ub,nonLinSolverOpts);
paramuq = uncertainty(nonlinfit);


% Residual vector function
% ------------------------------------------------------------------
% Function that provides vector of residuals, which is the objective
% function for the least-squares solvers
    function res = ResidualsFcn(p)
        
        % Non-linear model evaluation
        % ===============================
        A = Amodel(p);
        
        % Get condition of the operator
        if cond(A)>10
            illConditioned = true;
        else
            illConditioned = false;
        end
        
        % Regularization components
        % ===============================
        if illConditioned
            % Use an arbitrary axis
            ax = (1:1:size(A,2));
            % Get regularization operator
            RegOrder = min(size(A,2)-1,RegOrder);
            L = regoperator(ax,RegOrder);
            % If the parameter vector has not changed by much...
            if ~isempty(par_prev) && all(abs(par_prev-p)./p < alphaOptThreshold)
                % ...use the alpha optimized in the previous iteration
                alpha = regparam_prev;
            else
                % ...otherwise optimize with current settings
                alpha = selregparam(y,A,ax,'tikh',regparam,'RegOrder',RegOrder);
            end
            % Get components for LSQ fitting
            [AtAreg,Aty] = lsqcomponents(y,A,L,alpha,'tikhonov');
            
            % Store current iteration data for next one
            par_prev = p;
            regparam_prev = alpha;
        end
        
        
        if ~linearConstrained && ~illConditioned
            % Well-conditioned + Unconstrained
            % ====================================
            linfit = A\y;
            
            % Well-conditioned + Constrained
            % ====================================
        elseif linearConstrained && ~illConditioned
            linfit = linSolverFcn(A,y,lbl,ubl,linSolverOpts);
            
            
        elseif ~linearConstrained && illConditioned
            % Ill-conditioned + Unconstrained
            % ====================================
            linfit = AtAreg\Aty;
            
            
        elseif linearConstrained && illConditioned && ~nonNegativeOnly
            % Ill-conditioned + Constrained
            % ====================================
            linfit = linSolverFcn(AtAreg,Aty,lbl,ubl,linSolverOpts);
            
            
        elseif linearConstrained && illConditioned && nonNegativeOnly
            % Ill-conditioned + Non-Negativity
            % ====================================
            linfit = fnnls(AtAreg,Aty);
        end
        
        % Evaluate full model residual
        % ===============================
        res = A*linfit - y;
        
    end


% Uncertainty quantification
% ------------------------------------------------------------------
% Function that computes the covariance-based uncertainty quantification
% and returns the corresponding uncertainty structure
    function [paramuq] = uncertainty(parfit)
        
        
        if isempty(ubl)
            ubl = realmax*ones(numel(linfit),1);
        end
        if isempty(lbl)
            lbl = -realmax*ones(numel(linfit),1);
        end
        if isempty(ub)
            ub = realmax*ones(numel(linfit),1);
        end
        if isempty(lb)
            lb = -realmax*ones(numel(linfit),1);
        end
        
        % Compute the jacobian of the signal fit with respect to parameter set
        subidx_pnonlin = 1:numel(parfit);
        subidx_plin = numel(parfit)+[1:numel(linfit)];
        
        % Augmented Jacobian
        Jnonlin = jacobianest(@(p)Amodel(p)*linfit,parfit);
        Jlin = Amodel(parfit);
        if illConditioned
            Jreg = [zeros(size(L,1),numel(parfit)) regparam_prev*L];
        else
            Jreg = [];
        end
        J = [Jnonlin, Jlin; Jreg];
        
        % Suppress warnings for a moment
        warning('off','MATLAB:nearlySingularMatrix'), warning('off','MATLAB:singularMatrix')
        lastwarn('');
        
        % Estimate variance on experimental signal
        sigma2 = std(ResidualsFcn(parfit)).^2;
        
        % Estimate the covariance matrix by means of the inverse of Fisher information matrix
        covmatrix = sigma2.*inv(J.'*J);
        
        % Detect if there was a 'nearly singular' warning...
        [~, warnId] = lastwarn;
        if strcmp(warnId,'MATLAB:nearlySingularMatrix') || strcmp(warnId,'MATLAB:singularMatrix')
            % ...and if there was, then use a pseudoinverse instead of inverse
            covmatrix = sigma2.*sparse(pinv(full(J.'*J)));
            lastwarn('');
        end
        warning('on','MATLAB:nearlySingularMatrix'), warning('on','MATLAB:singularMatrix')
        
        % Construct uncertainty quantification structure for fitted parameters
        paramuq = uqst('covariance',[parfit(:); linfit(:)],covmatrix,[lb(:); lbl(:)],[ub(:); ubl(:)]);
    end


% Constrained linear least squares solver (toolbox-free)
% ------------------------------------------------------------------
% Alternative solver to lsqlin.m (Optimization toolbox). Solves the linear
% problem argmin_x ||Ax - y||^2 in its quadratic programming form, i.e.
% argmin_x 1/2 xt(2AtA)x+(-2Aty)tx using the open-source MINQ package:
% W. Huyer and A. Neumaier, MINQ8 (2017)
% General Definite and Bound Constrained Indefinite Quadratic Programming
    function [x,eval,exitflag] = lsqlin_QP(A,y,lbl,ubl,~)
        
        % MINQ is sensitive to vector orientation
        ubl = ubl(:);
        lbl = lbl(:);
        % Get Hessian
        H = 2*A.'*A;
        % Get linear term
        c = -2*A.'*y;
        
        % Unused settings
        gam = 0;
        print = 0;
        
        % Solve QP
        [x,eval,exitflag] = minq(gam,c,H,lbl,ubl,print);
    end
end



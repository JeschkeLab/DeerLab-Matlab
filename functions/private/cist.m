%
%   CIST Confidence interval structure constructor
%
%   cistruct = CIST('covariance',parfit,covmat,lb,ub)
%   Constructs the structure for covariance-based confidence intervals.
%
%   cistruct = CIST('bootstrap',samples)
%   Constructs the structure for bootstrapped confidence intervals.
%
%   Inputs:
%      parfit     N-element array of fitted values used as mean values in covariance-based CI
%      covmat     NxN-element covariance matrix
%      lb         N-element array of lower bounds used for fitting parfit
%      ub         N-element array of upper bounds used for fitting parfit
%      samples    MxN-element matrix of bootstrapped samples results

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function cistruct = cist(type,varargin)

%Parse inputs schemes
switch type
    
    case 'covariance'
        % Scheme 1: cist('covariance',parfit,covmat,lb,ub)
        parfit = varargin{1}(:);
        covmat = varargin{2};
        lb = varargin{3}(:);
        ub = varargin{4}(:);
        nParam = numel(parfit);
        
        if isempty(lb)
            lb = zeros(nParam,1) - realmax;
        end
        if isempty(ub)
            ub = zeros(nParam,1) + realmax;
        end
        
    case 'bootstrap'
        % Scheme 2: cist('bootstrap',samples)
        samples = varargin{1};
        nParam = numel(samples(1,:));
        
    otherwise
        error('Type not found. Must be: ''covariance'' or ''bootstrap''.')
end


%-------------------------------------------------------------------------
%Create confidence intervals structure
%-------------------------------------------------------------------------
switch type
    %Covariance-based CI specific fields
    case 'covariance'
        cistruct.mean = parfit;
        cistruct.median = parfit;
        cistruct.std = sqrt(diag(covmat));
        cistruct.covmat = covmat;
        cistruct.propagate = @(model,lb,ub)propagate(model,lb,ub);
        
        %Bootstrap-based CI specific fields
    case 'bootstrap'
        cistruct.mean = squeeze(mean(samples,1));
        cistruct.median = squeeze(median(samples,1));
        cistruct.std = squeeze(std(samples,[],1));
end

cistruct.percentile = @(p)percentile(p);
cistruct.ci = @(p)ci(p);
cistruct.pardist = @(n)pardist(n);

%-----------------------------------------------
% Parameter percentiles
%-----------------------------------------------
    function x = percentile(p)
        
        x = zeros(nParam,1);
        for n=1:nParam
            % Get parameter PDF
            dist = pardist(n);
            values = dist.values;
            pdf = dist.pdf;
            % Compute corresponding CDF
            cdf = cumsum(pdf);
            % Eliminate duplicates
            [cdf, index] = unique(cdf);
            % Interpolate requested percentile
            x(n) = interp1(cdf,values(index),p);
        end
    end

%-----------------------------------------------
% Covariance-based confidence intervals
%-----------------------------------------------
    function x = ci(coverage)
        
        alpha = 1 - coverage;
        p = 1 - alpha/2; % percentile
        
        switch type
            
            case 'covariance'
                % Compute covariance-based confidence intervals
                x(:,1) = max(lb,parfit - norm_inv(p)*sqrt(diag(covmat)));
                x(:,2) = min(ub,parfit + norm_inv(p)*sqrt(diag(covmat)));
                
            case 'bootstrap'
                p = (1 - coverage)/2;
                x(:,1) = percentile(p);
                x(:,2) = percentile(1-p);
        end
    end

%-----------------------------------------------
% Parameter distributions
%-----------------------------------------------
    function out = pardist(n)
        
        switch type
            
            case 'covariance'
                %Generate Gaussian distribution based on covariance matrix
                sig = sqrt(covmat(n,n));
                xmean = parfit(n);
                x = linspace(xmean-4*sig,xmean+4*sig,500);
                pdf = 1/sig/sqrt(2*pi)*exp(-((x-xmean)/sig).^2/2);
                
                %Clip the distributions at outside the boundaries
                pdf(x<lb(n)) = 0;
                pdf(x>ub(n)) = 0;
                
            case 'bootstrap'
                % Kernel density estimation
                [~,pdf,x] = kde(samples(:,n),500);
                
        end
        
        %Ensure normalization of the probability density function
        pdf = pdf/sum(pdf);
        
        %Store in structure
        out.pdf = pdf;
        out.values = x;
    end

%-----------------------------------------------
% Error Propagation (covariance-based only)
%-----------------------------------------------

    function modelcistruct = propagate(model,lb,ub)
        
        
        % Get jacobian of model to be propagated with respect to parameters
        jacobian = jacobianest(model,parfit);
        modelfit = model(parfit);
        
        % Clip at boundaries
        if isempty(lb)
            lb = zeros(numel(modelfit),1) - realmax;
        end
        if isempty(ub)
            ub = zeros(numel(modelfit),1) + realmax;
        end
        modelfit = max(modelfit,lb);
        modelfit = min(modelfit,ub);
        
        % Error progation
        modelcovmat = jacobian*covmat*jacobian.';
        
        % Construct new CI-structure for the model
        modelcistruct = cist(type,modelfit,modelcovmat,lb,ub);
        
    end


end
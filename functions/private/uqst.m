%
%   UQST Uncertainty quantification structure constructor
%
%   uqstruct = UQST('covariance',parfit,covmat,lb,ub)
%   Constructs the structure for covariance-based uncertainty quantificaton.
%
%   uqstruct = UQST('bootstrap',samples)
%   Constructs the structure for bootstrapped uncertainty quantificaton.
%
%   Inputs:
%      parfit     N-element array of fitted values used as mean values in covariance-based CI
%      covmat     NxN-element covariance matrix
%      lb         N-element array of lower bounds used for fitting parfit
%      ub         N-element array of upper bounds used for fitting parfit
%      samples    MxN-element matrix of bootstrapped samples results

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function uqstruct = uqst(type,varargin)

%Parse inputs schemes
switch type
    
    case 'covariance'
        % Scheme 1: uqst('covariance',parfit,covmat,lb,ub)
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
        % Scheme 2: uqst('bootstrap',samples)
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
        uqstruct.mean = parfit;
        uqstruct.median = parfit;
        uqstruct.std = sqrt(diag(covmat));
        uqstruct.covmat = covmat;
        
        %Bootstrap-based CI specific fields
    case 'bootstrap'
        means = squeeze(mean(samples,1));
        covmat = squeeze(samples).'*squeeze(samples)/size(samples,1) - means.*means.';
        uqstruct.mean = means;
        uqstruct.median = squeeze(median(samples,1));
        uqstruct.std = squeeze(std(samples,[],1));
        uqstruct.covmat = covmat;

end

uqstruct.percentile = @(p)percentile(p);
uqstruct.ci = @(p)ci(p);
uqstruct.pardist = @(n)pardist(n);
uqstruct.propagate = @(model,lb,ub)propagate(model,lb,ub);
uqstruct.type = type;

%-----------------------------------------------
% Parameter percentiles
%-----------------------------------------------
    function x = percentile(p)
  
        if p>100 || p<0
            error('The input must be a number between 0 and 100')
        end
        
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
            x(n) = interp1(cdf,values(index),p/100);
        end
    end

%-----------------------------------------------
% Covariance-based confidence intervals
%-----------------------------------------------
    function x = ci(coverage)
        
        if coverage>100 || coverage<0
            error('The input must be a number between 0 and 100')
        end
        
        alpha = 1 - coverage/100;
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
        
        if floor(n)-n~=0 || n>nParam || n<1
            error('The input must be an integer number between 1 and %i',nParam)
        end
        
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

    function modeluqstruct = propagate(model,lb,ub)
        
        if ~isa(model,'function_handle')
            error('The 1st input must be a valid function handle: @(parfit)model(__,parfit,__)')
        end
        
        % Evaluate model with fit parameters
        modelfit = model(parfit);

        % Validate input boundaries
        if isempty(lb)
            lb = zeros(numel(modelfit),1) - realmax;
        end
        if isempty(ub)
            ub = zeros(numel(modelfit),1) + realmax;
        end
        if numel(modelfit)~=numel(lb) || numel(modelfit)~=numel(ub)
            error('The 2nd and 3rd input arguments must have the same number of elements as the model output.')
        end
        
        % Get jacobian of model to be propagated with respect to parameters
        jacobian = jacobianest(model,parfit);
        
        % Clip at boundaries
        modelfit = max(modelfit,lb);
        modelfit = min(modelfit,ub);
        
        % Error progation
        modelcovmat = jacobian*covmat*jacobian.';
        
        % Construct new CI-structure for the model
        modeluqstruct = uqst(type,modelfit,modelcovmat,lb,ub);
        
    end


end
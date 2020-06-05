%
% FITMULTIMODEL Multi-model fitting of a distance distribution
%
%   P = FITMULTIMODEL(V,K,r,@basis,Nmax,method)
%   Fits a multi-model parametric distance distribution model to the dipolar
%   signal (V) using (@model) as the basis function, the dipolar kernel (K)
%   and distance axis (r). The function compares these multi-model distributions
%   with up to a maximum number of basis functions given by (Nmax) and 
%   determines the optimum one using the model selection criterion given in 
%   (method) ('AIC', 'BIC', or 'AICc'). The fitted distribution is returned in P.
%
%   P = FITMULTIMODEL(V,t,r,@basis,Nmax,method)
%   If a the default kernel is to be used, the time axis (t) can be passed
%   instead of the kernel.
%
%   P = FITMULTIMODEL({V1,V1,...},{K1,K2,K3},r,@basis,Nmax,method)
%   Passing multiple signals/kernels enables distance-domain global fitting
%   of the parametric models to single distributions. 
%   The multiple signals are passed as a cell array of arrays of sizes N1,N2,...
%   and a cell array of kernel matrices with sizes N1xM,N2xM,... must be 
%   passed as well.
%
%   P = FITMULTIMODEL({V1,V1,...},{t1,t2,t3},r,@basis,Nmax,method)
%   Similarly, time-domain global fitting can be used when passing time-domain
%   and the model time axes {t1,t2,...} of the corresponding signals.
%
%   [P,param,Pci,paramci,opt,metrics,Peval,stats] = FITMULTIMODEL(___)
%   If requested alongside the distribution (P), the optimal fit model
%   parameters (param), the optimal number of Gaussians (opt) and
%   evaluated selection metrics (metrics) are returned. The fitted distance
%   distributions fitted for ech multigauss model can be requested as a
%   fifth output argument (Peval). A structure containing different statistical
%   estimators of goodness of fit for the optimal model is returned as (stats).
%
%   P = FITMULTIMODEL(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
%   'Background' - Function handle to the background model to be fitted
%                 along the multigauss distance distribution model.
%                 Requires the time-axis to be passed instead of the
%                 kernel. Can be used with global fitting, where the same
%                 model will be applied to all signals.
%
%   'Lower' - Array [<r>_min FWHM_min] containing the lower bound for the
%             FWHM and mean distance of all the Gaussians.
%   'Upper' -  Array [<r>_max FWHM_max] containing the upper bound values
%              for the FWHM and mean distance of all the Gaussians.
%
%   See "help fitparamodel" for a detailed list of other property-value pairs
%   accepted by the function.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [Pfit,param,Pfitci,paramci,nGaussOpt,metrics,Peval,stats] = fitmultimodel(Vs,Ks,r,model,maxModels,method,varargin)


% Validate user input (S, K, r, and method are validated in lower-level functions)
if nargin<5
    error('Not enough input arguments.')
else
    validateattributes(maxModels,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'maxGaussians')
end
if nargin<6
    method = 'aicc';
elseif nargin==7
    varargin = [{method} varargin];
    method = 'aicc';
end

% Parse the optional parameters
%--------------------------------------------------------------
optionalProperties = {'Upper','Lower','Background','internal::parselater'};
[Upper,Lower,BckgModel] = parseoptional(optionalProperties,varargin);

% Control that the boundaries match the model and are appropiate
modelInfo = model();
nparam =  height(modelInfo);
paramNames = modelInfo.Parameter;
str = [];
for i=1:numel(paramNames), str = [str paramNames{i} ', ']; end, str(end-1:end) = ''; 
if ~isempty(Upper) && isempty(BckgModel) && length(Upper)~=nparam
    error('''Upper'' property must be a %i-element array of upper boundaries for the parameters [%s]',nparam,str)
elseif ~isempty(Upper) && ~isempty(BckgModel) && length(Upper)<(nparam+2)
    error('''Upper'' property must be a %i-element array [%s lambda_max Bparam_max]',nparam,str)
end
if ~isempty(Lower) && isempty(BckgModel)  && length(Lower)~=nparam
    error('''Lower'' property must be a %i-element array of lower boundaries for the parameters [%s]',nparam,str)
elseif ~isempty(Upper) && ~isempty(BckgModel) && length(Lower)<(nparam+2)
    error('''Lower'' property must be a %i-element array [%s lambda_min Bparam_min]',nparam,str)
end

% Remove used options from varargin so they are not passed to fitparamodel
for i = 1:numel(optionalProperties)
    Idx = find(cellfun(@(x)(ischar(x) && strcmpi(x,optionalProperties{i})),varargin));
    varargin(Idx:Idx+1) = [];
end

%Parse the required inputs for global fitting
if ~iscell(Vs)
   Vs = {Vs}; 
end
if ~iscell(Ks)
   Ks = {Ks}; 
end
for i=1:numel(Ks)
    if ~all(size(Ks{i}) > 1)
        ts{i} = Ks{i};
        Ks{i} = dipolarkernel(ts{i},r);
    end
end
if ~isempty(BckgModel) && ~exist('ts','var')
    error('Time axes must be provided for a time-domain fit.')
end
Nsignals = numel(Vs);

% Multi-component model construction
%--------------------------------------------------------------

% Compile list of multi-component models
multiModels = cell(maxModels,1);
if strcmp(func2str(model),'dd_gauss')
    % If basis function is a Gaussian, use built-in models
    multiModels{1} = @dd_gauss;
    if maxModels>=2, multiModels{2} = @dd_gauss2; end
    if maxModels>=3, multiModels{3} = @dd_gauss3; end
    if maxModels>=4, multiModels{4} = @dd_gauss4; end
    if maxModels>=5, multiModels{5} = @dd_gauss5; end
    for i = 6:maxModels
        multiModels{i} =  mixmodels(multiModels{i-1},@dd_gauss);
    end
elseif strcmp(func2str(model),'dd_rice')
    % If basis function is a Rician, use built-in models
    multiModels{1} = @dd_rice;
    if maxModels>=2, multiModels{2} = @dd_rice2; end
    if maxModels>=3, multiModels{3} = @dd_rice3; end
    if maxModels>=4, multiModels{4} = @dd_rice4; end
    if maxModels>=5, multiModels{5} = @dd_rice5; end
    for i = 6:maxModels
        multiModels{i} =  mixmodels(multiModels{i-1},@dd_rice);
    end
else
    % Otherwise mix the models
    multiModels{1} = model;
    for i = 2:maxModels
        multiModels{i} =  mixmodels(multiModels{i-1},model);
    end
end

% Preparation of parameter boundaries and start values
%--------------------------------------------------------------

% If the user has specified som boundaries then set the models boundaries appropiately
if ~isempty(Upper) || ~isempty(Lower)
    % Run over all multi-Gauss models
    for i = 1:maxModels
        % Get the info about the models
        info = multiModels{i}();
        modelnparam = height(info);
        boundary = zeros(1,modelnparam);
        paramNames = info.Parameter;

        %Get the indices of the different parameters on the mixed models
        paramidx = (1:nparam+1:modelnparam) + (0:nparam).';
        ampidx = paramidx(end,1:end-1);
        if ~isempty(Upper)
            for j=1:nparam
                boundary(paramidx(j,:)) = Upper(j);
            end
            boundary(ampidx) = 1; % Amplitudes upper bound
            if numel(Upper)>nparam+1
                boundary(numel(boundary)+1:numel(boundary) + numel(Upper(3:end))) = Upper(3:end); % Background upper bound
            end
            UpperBounds{i} = boundary;
        else
            UpperBounds = [];
        end
        boundary = zeros(1,modelnparam);
        if ~isempty(Lower)
            for j=1:nparam
                boundary(paramidx(j,:)) = Lower(j);
            end
            boundary(ampidx) = 0; % Amplitudes lower bound
            if numel(Lower)>nparam+1
                boundary(numel(boundary)+1:numel(boundary) + numel(Lower(3:end))) = Lower(3:end); %Background upper bound
            end
            LowerBounds{i} = boundary;
        else
            LowerBounds = [];
        end
               
    end
else
    % Otherwise just pass them empty to use the model defaults
    LowerBounds = [];
    UpperBounds = [];
end


% Preparation of the background model
%--------------------------------------------------------------
if ~isempty(BckgModel)
    if isempty(LowerBounds)
        for i = 1:maxModels
            info = multiModels{i}();
            range = [info.parameters(:).range];
            Plower = range(1:2:end-1);
            infoB = BckgModel();
            range = [infoB.parameters(:).range];
            Blower = range(1:2:end-1);
            LowerBounds{i} = [Plower repmat([0 Blower],1,Nsignals)];
            
        end
    end
    if isempty(UpperBounds)
        for i = 1:maxModels
            info = multiModels{i}();
            range = [info.parameters(:).range];
            Pupper = range(2:2:end);
            infoB = BckgModel();
            range = [infoB.parameters(:).range];
            Bupper = range(2:2:end);
            UpperBounds{i} = [Pupper repmat([1 Bupper],1,Nsignals)];
        end
    end
    for i = 1:maxModels
        DistModel = multiModels{i};
        info = DistModel();
        Nparam = height(info);
        Pparam = info.Start;
        infoB = BckgModel();
        Bparam = infoB.Start;
        lampars = Nparam + (1+numel(Bparam))*(1:Nsignals)-numel(Bparam);
        Bpars = lampars + 1;
        lam0 = 0.25;
        timeMultiGaussModels{i} = @(t,param,idx) (1 - param(lampars(idx)) + param(lampars(idx))*dipolarkernel(t,r)*DistModel(r,param(1:Nparam)) ).*BckgModel(t,param(Bpars(idx):Bpars(idx)+numel(Bparam)-1));
        param0{i} = [Pparam; repmat([lam0; Bparam],Nsignals,1)];
    end
end

% Optimal multi-component model selection
%--------------------------------------------------------------

% Run fitting and model selection to see which multi-Gauss model is optimal
if ~isempty(BckgModel)
    [nGaussOpt,metrics,fitparams,paramcis,stats] = selectmodel(timeMultiGaussModels,Vs,ts,method,param0,LowerBounds,UpperBounds,varargin);
else
    [nGaussOpt,metrics,fitparams,paramcis,stats] = selectmodel(multiModels,Vs,r,Ks,method,[],LowerBounds,UpperBounds,varargin);
end

% Calculate the distance distribution for the optimal multi-Gauss model
param = fitparams{nGaussOpt};
paramci = paramcis{nGaussOpt};
optModel = multiModels{nGaussOpt};
info = optModel();
nparam = height(info);
Pfit = optModel(r,param(1:nparam));
stats = stats{nGaussOpt};

% Uncertainty estimation
%--------------------------------------------------------------
if nargin>3
    %Loop over different signals
    lb = zeros(numel(r),1);
    Pfitci = paramci.propagate(@(par)optModel(r,par(1:nparam)),lb,[]);
    
%     jacobian = jacobianest(@(par)optModel(r,par),param(1:info.nparam));
%     % Calculate the confidence bands for the distance distribution
%     modelvariance = arrayfun(@(idx)full(jacobian(idx,:))*covmatrix*full(jacobian(idx,:)).',1:numel(r)).';
%     
%     Pfitci = cell(numel(critical),1);
%     for j=1:numel(critical)
%         upperci = Pfit + critical(j)*sqrt(modelvariance);
%         lowerci = max(0,Pfit - critical(j)*sqrt(modelvariance));
%         Pfitci{j} = [upperci(:) lowerci(:)];
%     end
%     %Do not return a cell if only one confidence level is requested
%     if numel(critical)==1
%         Pfitci = Pfitci{1};
%     end
end

if nargout>6
    Peval = zeros(maxModels,numel(r));
    for i = 1:maxModels
        info = multiModels{i}();
        nparam = height(info);
        p = fitparams{i};
        Peval(i,:) = multiModels{i}(r,p(1:nparam));
    end
end

return

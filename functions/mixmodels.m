%
% MIXMODELS Combine parametric models into one
%
%   newmodel = MIXMODELS(@model1,@model2,...,@modelN)
%   Combines the parametric model function handles (@model1,...,@modelN)
%   into a new parametric model function handle (newmodel). The models
%   must be passed as a cell array of function handles.
%
%   The parametric model function must be of the type as the models distributed
%   in DeerLab2. The returned function handle can be used for
%   parametric model fitting as the other models.
%
%   Example: dd_mix = MIXMODELS(@dd_gauss,@dd_gauss)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function mixModelFcn = mixmodels(varargin)

if numel(varargin)==1
    models = {varargin};
else
    models = varargin;
end

if numel(varargin)==0
    error('At least one model must be provided.')
end

if ~all(cellfun(@(M)isa(M,'function_handle'),models))
    error('Input arguments must all be function handles.')
end

% Detemine number of models to be mixed
nModels = numel(models);

% Combine the information structures of the models
%-------------------------------------------------------------------------------
% Add amplitudes for each model except last
for j = 1:nModels-1
    Info(j,1).Index = j;
    Info(j,1).Parameter = sprintf('Model %i: Amplitude',j);
    Info(j,1).Units = '  ';
    Info(j,1).Lower = 0;
    Info(j,1).Upper = 1;
    Info(j,1).Start = 1/nModels;
end
pidx_amp = 1:nModels-1;

% Combine info structures from all models
idx = pidx_amp(end);
for i = 1:nModels
    info = table2struct(models{i}());
    nparam = numel(info);
    pidx{i} = idx + (1:nparam);
    idx = idx + nparam;
    for j = 1:nparam
        info(j).Parameter = sprintf('Model %i: %s',i,info(j).Parameter);
    end
    Info = [Info; info];
end

% Mixed model function handle
%-------------------------------------------------------------------------------
mixModelFcn = @mixedFunction;

    % Function to allow request of information structure or model values
  function output = mixedFunction(varargin)

        if nargin==0
            output = struct2table(Info);
            return
        end
        
        if nargin<2
            error('At least two input arguments (r,param) are required.')
        elseif  nargin>3
            error('Only two input arguments are allows.')
        end
        
        x = varargin{1};
        params = varargin{2};
        if ~iscolumn(x)
            x = x.';
        end
        
        amp = params(pidx_amp);
        amp(end+1) = max(1-sum(amp),0);
        
        y = 0;
        for k = 1:numel(models)
            y = y + amp(k)*models{k}(x,params(pidx{k}));
        end
        output = y;
        
    end
    
end

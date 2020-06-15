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
% Combine info structures from all models
idx = 0;
Info = [];
for i = 1:nModels
    info = models{i}();
    nparam = numel(info);
    pidx{i} = idx + (1:nparam);
    idx = idx + nparam;
    for j = 1:nparam
        info(j).Index =  pidx{i}(j);
        info(j).Parameter = sprintf('Model %i: %s',i,info(j).Parameter);
    end
    
    % Add amplitudes for each model
    ampInfo.Index = idx+1;
    ampInfo.Parameter = sprintf('Model %i: Amplitude',i);
    ampInfo.Units = '  ';
    ampInfo.Lower = 0;
    ampInfo.Upper = 1;
    ampInfo.Start = 1/nModels;
    Info = [Info info ampInfo];
    pidx_amp(i) = numel(Info);
    idx = idx + 1;
end

% Mixed model function handle
%-------------------------------------------------------------------------------
mixModelFcn = @mixedFunction;

    % Function to allow request of information structure or model values
  function output = mixedFunction(varargin)

        if nargin==0
            output = Info;
            return
        end
        
        if nargin<2
            error('At least two input arguments (axis,parameters) are required.')
        elseif  nargin>3
            error('Only two input arguments are allowed.')
        end
        
        x = varargin{1};
        params = varargin{2};
        if ~iscolumn(x)
            x = x.';
        end
        
        amp = params(pidx_amp);        
        y = 0;
        for k = 1:numel(models)
            y = y + amp(k)*models{k}(x,params(pidx{k}));
        end
        output = y;
        
    end
    
end

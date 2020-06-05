%
% PARAMODEL Converts a function handle to a DeerLab parametric model
%
%       g = PARAMODEL(f,n)
%       Converts the input function handle (f) to a valid DeerLab
%       parametric model compatible with the fit functions. The number of
%       parameters in the resulting parametric model must be specified by
%       (n). The initial guess values of all parameters are set to zero.
%       The resulting parametric model g is returned in as a function
%       handle.
%
%       g = PARAMODEL(f,param0)
%       If an array of inital guess values are passed, these are set into
%       the resulting parametric model. The number of parameters is computed
%       from the length of the array.
%
%       g = PARAMODEL(f,param0,lower,upper)
%       Lower and upper bounds for the parameters can be provided as optional
%       input arguments.

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function model = paramodel(handle,param0,lower,upper,normalize)

if nargin<2
  error('At least two inputs (function handle and parameter vector) are needed.');
end
if nargin<4
   normalize = false; 
end

nParam = numel(param0);

if nargin<3 || isempty(lower)
    lower = zeros(nParam,1) - realmax;
    upper = zeros(nParam,1) + realmax;
end
if nargin<4 || isempty(upper)
    upper = zeros(nParam,1) + realmax;
end

if length(param0)~=length(upper) || length(param0)~=length(lower)
   error('The inital guess and upper/lower boundaries must have equal length.') 
end

model = @myparametricmodel;

% Define the raw structure of the DeerLab parametric model functions
    function Output = myparametricmodel(ax,param,idx)
        
        if nargin==0
            % If no inputs given, return info about the parametric model
            for i = 1:nParam
                info(i).Index = i;
                info(i).Parameter = ' ';
                info(i).Units = ' ';
                info(i).Lower = lower(i);
                info(i).Upper = upper(i);
                info(i).Start = param0(i);
            end
            
            Output = struct2table(info);
            
        elseif nargin >= 2
            
            % If user passes them, check that the number of parameters matches the model
            if length(param)~=nParam
                error('The number of input parameters does not match the number of model parameters.')
            end
            
            % If necessary inputs given, compute the model distance distribution
            if nargin(handle)==3
                Output = handle(ax,param,idx);
            else
                Output = handle(ax,param);
            end
            if normalize
                Output = Output/sum(Output)/mean(diff(ax));
            end
            if ~iscolumn(Output)
                Output = Output.';
            end
        else
            % Else, the user has given wrong number of inputs
            error('Model requires two input arguments.')
        end
        
    end

end

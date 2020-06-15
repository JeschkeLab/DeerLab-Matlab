%
% BOOTAN Bootstrap analysis for uncertainty estimation
%
%   [bootci,stats] = BOOTAN(fcn,V,Vfit)
%   [bootci,stats] = BOOTAN(fcn,V,Vfit,Nsamples)
%   [bootci,stats] = BOOTAN(fcn,{V1,V2,___},{Vfit1,Vfit2,___},Nsamples)
%   [bootci,stats] = BOOTAN(___,'Name',value,___)
%
%   Performs a bootstrap analysis of the output variables returned by the
%   function (fcn) given an experimental signal (V) and its fit (Vfit).
%
%   Inputs:
%     fcn        function handle to the data analysis function, must be
%                accept V as sole input: fcn = @(V) myfunction(___,V,___)
%     V          experimental signal
%     Vfit       fitted signal
%     Nsamples   number of bootstrap samples to evaluate (default=1000)
%
%  Name-value pairs:
%
%     'Verbose'    - Display progress information on the command window (true/false)
%     'Resampling' - Resampling method for bootstrap samples
%                        'gaussian' - sample noise from a Gaussian distribution (default)
%                        'residual' - sample noise from the fit residuals
%
%   Outputs:
%     bootci    structure array with confidence intervals for each variable
%       .ci99   99%-Confidence intervals (percentile-based)
%       .ci95   95%-Confidence intervals (percentile-based)
%       .ci50   50%-Confidence intervals (percentile-based)
%     stats     structure array with summary statistics for each variable
%       .median medians of the output variables
%       .mean   means of the output variables
%       .std    standard deviations of the output variables
%       .p1     2nd  percentiles of the output variables
%       .p25    25th percentiles of the output variables
%       .p75    75th percentiles of the output variables
%       .p99    98th percentiles of the output variables
%       .boothist    structure containing histogram data of the variable distributions
%                        .edges  Histogram edges (X-data)
%                        .bins   Histogram bins (in probability) (Y-data)
%       .bootdist    structure containing the kernel-density estimation
%                    distribution of the histogram
%                        .values   Evaluated variable values
%                        .pdf      Estimated probability density function
%
%   The .bootdist and .boothist fields are only computed for non-vectorial
%   variables, e.g. model parameters, but not for vectorial variables, e.g.
%   non-parametric distance distributions, backgrounds,...
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function stats = bootan(fcn,Vexp,Vfit,nSamples,varargin)

% Parse input scheme: bootan(fcn,Vexp,Vfit,'prop',value,___)
if nargin>=4 && ischar(nSamples)
    varargin = [{nSamples} varargin];
    nSamples = [];
end

% Parse the optional parameters in varargin
optionalProperties = {'Verbose','Resampling'};
[Verbose,Resampling] = parseoptions(varargin);

if nargin<4 || isempty(nSamples)
    nSamples = 1000;
else
    validateattributes(nSamples,{'numeric'},{'nonnegative','scalar','nonempty'},mfilename,'Nsamples')
end
if isempty(Resampling)
    Resampling = 'gaussian';
else
    validateattributes(Resampling,{'char'},{'nonempty'},mfilename,'Resampling')
    validInputs = {'residual','gaussian'};
    Resampling = validatestring(Resampling,validInputs);
end

if isempty(Verbose)
    Verbose = false;
else
    validateattributes(Verbose,{'logical'},{'nonempty'},mfilename,'Verbose')
end

if ~iscell(Vexp)
    Vexp = {Vexp};
end
if ~iscell(Vfit)
    Vfit = {Vfit};
end
if numel(Vexp) ~= numel(Vfit)
    error('The same number of signals V and fits Vfit must be provided.')
end
nSignals = numel(Vexp);
for i=1:numel(Vfit)
    if numel(Vexp{i}) ~= numel(Vfit{i})
        error('V{%i} and Vfit{%i} must have the same number of element.',i,i)
    end
end
if ~isa(fcn,'function_handle')
    error('The 1st input must be a valid function handle of the form @(V)fcn(__,V,__)')
end

% Determine the number of outputs returned by function handle and allocate arrays
nout = getmaxnargout(fcn,Vexp);
evals = cell(1,nout);
sizeOut = cell(1,nout);

% Get residuals and estimate standard deviation
residuals = cell(1,nSignals);
sigma = zeros(1,nSignals);
for i=1:nSignals
    residuals{i} = Vfit{i} - Vexp{i};
    sigma(i) = std(residuals{i});
end
% Generate bootstrap samples from data
%-------------------------------------------------------------------------------
Vresample = cell(1,nSignals);
warning('off','all')
for iSample = 1:nSamples
    
    % Inform of progress if requested
    if Verbose
        S = sprintf('Bootstrapping: %i/%i samples finished',iSample,nSamples);
        fprintf(S);
    end
    
    for i=1:nSignals
        % Get a bootstrap sample
        if iSample>1
            
            %Determine requested re-sampling method
            switch Resampling
                case 'gaussian'
                    % Resample from a Gaussian distribution with
                    % variance estimated from the residuals
                    Vresample{i} = Vfit{i} + sigma(i)*randn(size(Vfit{i}));
                case 'residual'
                    % Resample from the residual directly
                    Vresample{i} =  Vfit{i} + residuals{i}(randperm(numel(Vfit{i})));
            end
        else
            % Use original data on the first run
            Vresample{i} = Vexp{i};
        end
    end
    % Run the model function with bootstrap sample
    varargsout = cell(1,nout);
    [varargsout{:}] = fcn(Vresample);
    
    % Assert that all outputs are strictly numerical
    numericOutput = cellfun(@(x)isnumeric(x),varargsout);
    if ~all(numericOutput)
        error('Non-numeric output arguments by the input function are not accepted.');
    end
    
    % Loop over the different outputs of the function
    for iOut = 1:nout
        out = varargsout{iOut};
        
        if isempty(sizeOut{iOut})
            sizeOut{iOut} = size(out);
        elseif any(sizeOut{iOut}~=size(out))
            error(['Inconsistent output variable size. ',...
                'One of the outputs of the analyzed function is changing its size in between runs. ',...
                'To solve this, fix the axis of the output and interpolate the result. \n%s'],...
                '  Ex/ outFix = interp1(varAxis,out,fixAxis,''pchip'')')
        end
        % Convert vectors to rows
        if iscolumn(out)
            out = out.';
        end
        nParam(iOut) = numel(out);
        % Store outputs in an N-dimensional array
        evals{iOut} = cat(1,evals{iOut},shiftdim(out,-1));
    end
    
    %Reset the printed line (for one-line updat printing)
    if Verbose
        fprintf(repmat('\b',1,numel(S)));
    end
end
warning('on','all')

% Compile statistics for all parameters from bootstrap samples
%-------------------------------------------------------------------------------
for iOut = 1:numel(evals)
    boots = evals{iOut};
    stats{iOut} = uqst('bootstrap',boots);
end

end

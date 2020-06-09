%
% EX_4PDEER Single-pathway 4-pulse DEER experiment model 
%
%   info = EX_4PDEER()
%   Returns an (info) table of model parameters and boundaries.
%
%   pathways = EX_4PDEER(param)
%   Computes the dipolar pathway information array according to the paramater
%   array (param).
%
%
% PARAMETERS
%    -----------------------------------------------------------------
%     Index  Parameter           Units    Lower    Upper    Start
%    -----------------------------------------------------------------
%       1    Modulation depth               0        1       0.3 
%    -----------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = ex_4pdeer(param)

nParam = 1;

if nargin>1
    error('Model requires one input argument.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = "Modulation depth";
    info(1).Units = '  ';
    info(1).Lower = 0;
    info(1).Upper = 1;
    info(1).Start = 0.3;
    
    output = struct2table(info);
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Dipolar pathways
lambda = param(1);
pathway(1,:) = [1-lambda NaN];
pathway(2,:) = [lambda 0];
output = pathway;

end

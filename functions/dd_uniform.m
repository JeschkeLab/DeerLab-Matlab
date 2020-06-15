%
% DD_UNIFORM Uniform distribution parametric model
%
%   info = DD_UNIFORM
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_UNIFORM(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    --------------------------------------------------
%     Index  Parameter   Units   Lower   Upper   Start
%    --------------------------------------------------
%       1    Left edge    nm      0.1     6       2.5 
%       2    Right edge   nm      0.2     20       3 
%    --------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_uniform(r,param)

nParam = 2;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Left edge';
    info(1).Units = 'nm';
    info(1).Lower = 0.1;
    info(1).Upper = 6;
    info(1).Start = 2.5;
    
    info(2).Index = 2;
    info(2).Parameter = 'Right edge';
    info(2).Units = 'nm';
    info(2).Lower = 0.2;
    info(2).Upper = 20;
    info(2).Start = 3;
    
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
  error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
rL = min(abs(param));
rR = max(abs(param));
P = zeros(numel(r),1);
P(r>=rL & r<=rR) = 1;
P = P/trapz(r,P);

output = P;

return
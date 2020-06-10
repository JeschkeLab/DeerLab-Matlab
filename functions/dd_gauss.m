%
% DD_GAUSS Gaussian distribution parametric model
%
%   info = DD_GAUSS
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_GAUSS(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    -----------------------------------------------------------
%     Index    Parameter     Units     Lower    Upper    Start
%    -----------------------------------------------------------
%       1       Center        nm        1        20       3.5 
%       2       FWHM          nm       0.2       5        0.5 
%    -----------------------------------------------------------

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_gauss(r,param)

nParam = 2;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Center';
    info(1).Units = 'nm';
    info(1).Lower = 1;
    info(1).Upper = 20;
    info(1).Start = 3.5;
    
    info(2).Index = 2;
    info(2).Parameter = 'FWHM';
    info(2).Units = 'nm';
    info(2).Lower = 0.2;
    info(2).Upper = 5;
    info(2).Start = 0.5;
    
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
fwhm = param(2);
r0 = param(1);
a = 1;
P = multigaussfun(r,r0,fwhm,a);

output = P;

return
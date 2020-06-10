%
% DD_GENGAUSS Generalized Gaussian distribution parametric model
%
%   info = DD_GENGAUSS
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_GENGAUSS(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    ---------------------------------------------------
%     Index  Parameter  Units  Lower   Upper   Start
%    ---------------------------------------------------
%       1    Center     nm       1      20      3.5 
%       2    FWHM       nm      0.2      5      0.5 
%       3    Kurtosis           0.25    15       5 
%    ---------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_gengauss(r,param)

nParam = 3;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %I f no inputs given, return info about the parametric model
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
    
    info(3).Index = 3;
    info(3).Parameter = 'Kurtosis';
    info(3).Units = '  ';
    info(3).Lower = 0.25;
    info(3).Upper = 15;
    info(3).Start = 5;
    
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
r0 = param(1);
width = param(2);
beta = param(3);

sigma = width/(2*sqrt(2*log(2)));
x = abs(r(:)-r0)/sigma;
P = beta/(2*sigma*gamma(1/beta))*exp(-x.^beta);

% Normalize
if ~all(P==0)
    P = P/trapz(r,P);    
end

output = P;

return
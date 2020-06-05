
% DD_RICE 3D-Rice distribution parametric model
%
%   info = DD_RICE
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_RICE(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
%    ----------------------------------------------------
%     Index  Parameter   Units   Lower    Upper    Start
%    ----------------------------------------------------
%       1    Center       nm       1       10       3.5 
%       2    Width        nm      0.1       5       0.7 
%    ----------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = dd_rice(r,param)

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
    info(1).Upper = 10;
    info(1).Start = 3.5;
    
    info(2).Index = 2;
    info(2).Parameter = 'Width';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 5;
    info(2).Start = 0.7;
    
    output = struct2table(info);
    return
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute non-central chi distribution with 3 degrees of freedom (a 3D Rician)
nu = param(1);
sig = param(2);
a = 1;
P = multirice3d(r,nu,sig,a);

output = P;

return

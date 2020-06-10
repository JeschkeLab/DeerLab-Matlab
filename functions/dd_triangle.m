%
% DD_TRIANGLE Triangle distribution parametric model
%
%   info = DD_TRIANGLE
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_TRIANGLE(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    ----------------------------------------------------
%     Index  Parameter   Units    Lower   Upper    Start
%    ----------------------------------------------------
%       1    Center       nm       1       20       3.5 
%       2    Width left   nm      0.1       5       0.3 
%       3    Width right  nm      0.1       5       0.3 
%    ----------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_triangle(r,param)

nParam = 3;

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
    info(2).Parameter = 'Width left';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 5;
    info(2).Start = 0.3;
    
    info(3).Index = 3;
    info(3).Parameter = 'Width right';
    info(3).Units = 'nm';
    info(3).Lower = 0.1;
    info(3).Upper = 5;
    info(3).Start = 0.3;
    
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
wL = abs(param(2));
wR = abs(param(3));
rL = r0 - wL;
rR = r0 + wR;
idxL = r>=r0-wL & r<=r0;
idxR = r<=r0+wR & r>=r0;
P = zeros(numel(r),1);
if wL>0
    P(idxL) = (r(idxL)-rL)/wL/(wL+wR);
end
if wR>0
    P(idxR) = -(r(idxR)-rR)/wR/(wL+wR);
end

if any(P~=0)
    P = P/trapz(r,P);
end

output = P;

return
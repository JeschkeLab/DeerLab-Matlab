%
% BG_EXP Exponential background model
%
%   info = BG_EXP
%   Returns an (info) table of model parameters and boundaries.
%
%   B = BG_EXP(t,param)
%   B = BG_EXP(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
%    ---------------------------------------------------
%     Index  Parameter    Units   Lower   Upper   Start
%    ---------------------------------------------------
%       1    Decay Rate   us^-1     0      200    0.35 
%    ---------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.



function output = bg_exp(t,param,lambda)

nParam = 1;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model  
    info(1).Index = 1;
    info(1).Parameter = 'Decay Rate';
    info(1).Units = 'us^-1';
    info(1).Lower = 0;
    info(1).Upper = 200;
    info(1).Start = 0.35;
    
    output = info;
    return
end

if nargin<3
    lambda = 1;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% If necessary inputs given, compute the model background
kappa = param(1);
B = exp(-lambda*kappa*abs(t));
B = B(:);
output = B;


return
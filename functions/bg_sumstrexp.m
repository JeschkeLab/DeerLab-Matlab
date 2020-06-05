%
% BG_SUMSTREXP Sum of two stretched exponentials background model
%
%   info = BG_SUMSTREXP
%   Returns an (info) table of model parameters and boundaries.
%
%   B = BG_SUMSTREXP(t,param)
%   B = BG_SUMSTREXP(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
%    ---------------------------------------------------------------------------
%     Index  Parameter                        Units    Lower    Upper    Start
%    ---------------------------------------------------------------------------
%       1    Decay Rate of 1st component      us^-1      0       200      0.25 
%       2    Stretch Factor of 1st component             0         6       1 
%       3    Amplitude of 1st component                  0         1      0.5 
%       4    Decay Rate of 2nd component      us^-1      0       200      0.25 
%       5    Stretch Factor of 2nd component             0         6       1 
%    ---------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = bg_sumstrexp(t,param,lambda)

nParam = 5;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Decay Rate of 1st component';
    info(1).Units = 'us^-1';
    info(1).Lower = 0;
    info(1).Upper = 200;
    info(1).Start = 0.25;
    
    info(2).Index = 2;
    info(2).Parameter = 'Stretch Factor of 1st component';
    info(2).Units = '     ';
    info(2).Lower = 0;
    info(2).Upper = 6;
    info(2).Start = 1;
  
    info(3).Index = 3;
    info(3).Parameter = 'Amplitude of 1st component';
    info(3).Units = '     ';
    info(3).Lower = 0;
    info(3).Upper = 1;
    info(3).Start = 0.5;
    
    info(4).Index = 4;
    info(4).Parameter = 'Decay Rate of 2nd component';
    info(4).Units = 'us^-1';
    info(4).Lower = 0;
    info(4).Upper = 200;
    info(4).Start = 0.25;
    
    info(5).Index = 5;
    info(5).Parameter = 'Stretch Factor of 2nd component';
    info(5).Units = '     ';
    info(5).Lower = 0;
    info(5).Upper = 6;
    info(5).Start = 1;
    
    output = struct2table(info);
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
kappa1 = param(1);
d1 = param(2);
w1 = param(3);
kappa2 = param(4);
d2 = param(5);
strexp1 = exp(-lambda*kappa1*abs(t).^d1);
strexp2 = exp(-lambda*kappa2*abs(t).^d2);
B = w1*strexp1 + (1-w1)*strexp2;
B = B(:);
output = B;


return
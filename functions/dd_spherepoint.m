%
% DD_SPHEREPOINT One particle distanced from particles distributed on a
% sphere 
%
%   info = DD_SPHEREPOINT
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_SPHEREPOINT(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
%    --------------------------------------------------------
%     Index  Parameter        Units   Lower   Upper   Start
%    -------------------------------------------------------
%       1    Sphere radius      nm     0.1     20      1.5 
%       1    Distance to point  nm     0.1     20      3.5 
%    --------------------------------------------------------
%
%   See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
%        http://doi.org/10.1016/j.jmr.2013.01.007
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_spherepoint(r,param)

nParam = 2;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0    
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Sphere radius';
    info(1).Units = 'nm';
    info(1).Lower = 0.1;
    info(1).Upper = 20;
    info(1).Start = 1.5;
  
    info(2).Index = 2;
    info(2).Parameter = 'Distance to point';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 20;
    info(2).Start = 3.5;
    
    output = struct2table(info);
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
  error('The number of input parameters does not match the number of model parameters.')
end

%Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
R = param(1);
d = param(2);
P = zeros(numel(r),1);
idx = r >= d - R & r<= d + R; 
P(idx) = 3*r(idx).*(R^2 - (d - r(idx)).^2)./(4*d*R.^3);

if ~all(P==0)
P = P/trapz(r,P);
end

output = P;

return
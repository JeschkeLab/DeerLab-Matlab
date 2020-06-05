%
% DD_SHELL Uniform spherical shell
%
%   info = DD_SHELL
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_SHELL(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    --------------------------------------------------------
%     Index  Parameter       Units  Lower    Upper    Start
%    --------------------------------------------------------
%       1    Shell radius    nm     0.1      20       1.5 
%       2    Shell thickness  nm     0.1      20       0.5 
%    --------------------------------------------------------
%
%   See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
%        http://doi.org/10.1016/j.jmr.2013.01.007
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_shell(r,param)

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
    info(2).Parameter = 'Shell thickness';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 20;
    info(2).Start = 0.5;
    
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
R1 = param(1);
w = param(2);
R2 = R1 + w;
P = zeros(numel(r),1);

P = R2^6*pb(r,R2) - R1^6*pb(r,R1) - 2*(R2^3 - R1^3)*pbs(r,R1,R2);
P = P/(R2^3 - R1^3)^2;

if ~all(P==0)
P = P/trapz(r,P);
end

output = P;


function P = pb(r,R) 

idx = r >= 0 & r<= 2*R;
P = zeros(numel(r),1);
P(idx) = 3*r(idx).^5/(16*R^6) - 9*r(idx).^3/(4*R^4) + 3*r(idx).^2/(R^3);

end

function P = pbs(r,R1,R2)

P = zeros(numel(r),1);
%Case1
idx = r >= 0 & r < min(2*R1,R2 - R1); 
P(idx) = 12*r(idx).^3*R1^2 - r(idx).^5;

%Case2
idx = r >= R2 - R1 & r < 2*R1;
P(idx) = 8*r(idx).^2*(R2^3 - R1^3) - 3*r(idx)*(R2^2 - R1^2)^2 - 6*r(idx).^3*(R2 - R1)*(R2 + R1);

%Case3
idx = r >= 2*R1 & r < R2 - R1;
P(idx) = 16*r(idx).^2*R1^3;

%Case4
idx = r >= max(R2 - R1,2*R1) & r < R1 + R2;
P(idx) = r(idx).^5 - 6*r(idx).^3*(R2^2 + R1^2) + 8*r(idx).^2*(R2^3 + R1^3) - 3*r(idx)*(R2^2 - R1^2)^2;

P = P*3/(16*R1^3*(R2^3 - R1^3));

end

end

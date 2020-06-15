%
% DD_SHELLVOIDSHELL Particles distributed on a sphere inside a spherical shell separated by a void 
%
%   info = DD_SHELLVOIDSHELL
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_SHELLVOIDSHELL(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    -------------------------------------------------------------
%     Index  Parameter              Units   Lower   Upper   Start
%    -------------------------------------------------------------
%       1    Inner sphere radius}    nm      0.1     20     0.75 
%       2    Inner shell thickness   nm      0.1     20      1 
%       3    Outer shell thickness   nm      0.1     20      1 
%       4    Shell-shell separation  nm      0.      2      0.5 
%    -------------------------------------------------------------
%
%   See: D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63 
%        http://doi.org/10.1016/j.jmr.2013.01.007
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_shellvoidshell(r,param)

nParam = 4;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Inner sphere radius';
    info(1).Units = 'nm';
    info(1).Lower = 0.1;
    info(1).Upper = 20;
    info(1).Start = 0.75;
    
    info(2).Index = 2;
    info(2).Parameter = 'Inner shell thickness';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 20;
    info(2).Start = 1.0;
 
    info(3).Index = 3;
    info(3).Parameter = 'Outer shell thickness';
    info(3).Units = 'nm';
    info(3).Lower = 0.1;
    info(3).Upper = 20;
    info(3).Start = 1.0;
    
    info(4).Index = 4;
    info(4).Parameter = 'Shell-shell separation';
    info(4).Units = 'nm';
    info(4).Lower = 0.1;
    info(4).Upper = 20;
    info(4).Start = 0.5;
    
    output = info;
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
w1 = param(2);
w2 = param(3);
d = param(4);

R2 = R1 + w1;
R3 = R1 + w1 + w2;
R4 = R1 + w1 + w2 + d;

delta21 = R2^3 - R1^3;
delta31 = R3^3 - R1^3;
q31 = delta31*pbs(r,R1,R3);
delta32 = R3^3 - R2^3;
q32 = delta32*pbs(r,R2,R3);
delta41 = R4^3 - R1^3;
q41 = delta41*pbs(r,R1,R4);
delta42 = R4^3 - R2^3;
q42 = delta42*pbs(r,R2,R4);
delta43 = R4^3 - R3^3;

P = (R1^3*(q31 - q41) + R2^3*(q42 - q32))/(delta43*delta21);
P = round(P,15);

if ~all(P==0)
P = P/trapz(r,P);
end

output = P;

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
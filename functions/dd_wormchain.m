function output = dd_wormchain(r,param)
%
% WORMCHAIN Worm-like chain model near the rigid limit
%
%   info = WORMCHAIN
%   Returns an (info) table of model parameters and boundaries.
%
%   P = WORMCHAIN(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
%    ----------------------------------------------------------
%     Index  Parameter           Units   Lower   Upper   Start
%    ----------------------------------------------------------
%       1    Chain length         nm      1.5     20      3.7 
%       2    Persistence length   nm       2      100     10 
%    ----------------------------------------------------------
%
% See: J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
%      https://doi.org/10.1103/PhysRevLett.77.2581
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


nParam = 2;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Chain length';
    info(1).Units = 'nm';
    info(1).Lower = 1.5;
    info(1).Upper = 20;
    info(1).Start = 3.7;
    
    info(2).Index = 2;
    info(2).Parameter = 'Persistence length';
    info(2).Units = 'nm';
    info(2).Lower = 2;
    info(2).Upper = 100;
    info(2).Start = 10;
    
    output = info;
    return
end

%If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

%Parse input
if ~iscolumn(r)
    r = r.';
end
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

%Prepare parameters
L = param(1);
Lp = param(2);
kappa = Lp/L;


%Get normalized distance axis
normDistAxis=r/L;
P=zeros(size(r));
crit=kappa*(1 - normDistAxis);
%Compute ditribution using two terms of the expansion

rcrit = normDistAxis(crit>0.2);
P(crit>0.2)=2*kappa/(4*pi)*(pi^2*(-1)^(2)*exp(-kappa*pi^2*(1-rcrit)) ...
    + pi^2*4*(-1)^(3)*exp(-kappa*pi^2*4*(1-rcrit)));

rcrit = normDistAxis(crit>0 & crit<0.2);
P(crit>0 & crit<0.2) = kappa/(4*pi*2*sqrt(pi))*((1./(kappa*(1 - rcrit)).^(3/2).*exp(-(1 - 1/2)^2./(kappa*(1 - rcrit))).*(4.*((1 - 1/2)./sqrt(kappa*(1-rcrit))).^2-2)) ...
    + 1./(kappa*(1 - rcrit)).^(3/2).*exp(-(2 - 1/2)^2./(kappa*(1 - rcrit))).*(4.*((2 - 1/2)./sqrt(kappa*(1-rcrit))).^2-2));

%Normalize integral
if ~all(P==0)
P = P/trapz(r,P);
end
output = P;


return

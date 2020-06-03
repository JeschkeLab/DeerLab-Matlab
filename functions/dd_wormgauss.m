function output = dd_wormgauss(r,param)
%
% WORMGAUSS Worm-like chain model near the rigid limit with Gaussian
%           convolution
%
%   info = WORMGAUSS
%   Returns an (info) structure containing the specifics of the model.
%
%   P = WORMGAUSS(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  L      3.7     1.5         10             length of the worm-like chain
% param(2)  Lp     10      2           100            persistence length
% param(3)  sigma  0.2     0.001       2              Gaussian standard deviation 
% --------------------------------------------------------------------------
%
% See: J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
%      https://doi.org/10.1103/PhysRevLett.77.2581
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


nParam = 3;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Worm-like chain model near rigid limit';
    info.nparam  = nParam;
    info.parameters(1).name = 'Chain length';
    info.parameters(1).range = [1.5 10];
    info.parameters(1).default = 3.7;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'Persistence length';
    info.parameters(2).range = [2 100];
    info.parameters(2).default = 10;
    info.parameters(2).units = 'nm';
 
    info.parameters(3).name = 'Gaussian standard deviation';
    info.parameters(3).range = [0.001 2];
    info.parameters(3).default = 0.2;
    info.parameters(3).units = 'nm';
    
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
L=param(1);
Lp=param(2);
sigma = param(3);

kappa=Lp/L;


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

%Compute Gaussian convolution window
[~,idx] = max(P);
gauss = exp(-((r - r(idx))./sigma).^2);

%Convolution with size retention
P = conv(gauss,P,'full');

%Adjust new convoluted axis
rconv = linspace(min(r),max(r)*2,numel(P));
[~,idxconv] = max(P);
rconv = rconv  - abs(r(idx) - rconv(idxconv));

%Interpolate down to original axis
P = interp1(rconv,P,r,'pchip');

%Normalize integral
if ~all(P==0)
P = P/trapz(r,P);
end
output = P;

return

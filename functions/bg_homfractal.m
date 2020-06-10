%
% BG_HOMFRACTAL Multi-pulse DEER background in a fractal medium
%
%   info = BG_HOMFRACTAL
%   Returns an (info) table of model parameters and boundaries.
%
%   B = BG_HOMFRACTAL(t,param)
%   B = BG_HOMFRACTAL(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
%    --------------------------------------------------------------
%     Index  Parameter           Units      Lower    Upper   Start
%    --------------------------------------------------------------
%       1    Spin concentration  umol/dm^d   0.01    5000     50  
%       2    Dimension                        0      6        3  
%    --------------------------------------------------------------


% name    symbol default lower bound upper bound
% ------------------------------------------------------------------
% PARAM(1)  c       50        0.01       5000     fractal spin concentration (umol/dm^d)
% PARAM(2)  d       3          0         6        fractal dimension
% ------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = bg_homfractal(t,param,lambda)

nParam = 2;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Spin concentration';
    info(1).Units = 'umol/dm^d';
    info(1).Lower = 0.01;
    info(1).Upper = 5000;
    info(1).Start = 50;
    
    info(2).Index = 2;
    info(2).Parameter = 'Dimension';
    info(2).Units = '         ';
    info(2).Lower = eps;
    info(2).Upper = 6-6*eps;
    info(2).Start = 3;
    
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
conc = param(1); % umol/dm^d
d = param(2);

NA = 6.02214076e23; % Avogadro constant, mol^-1
conc = conc*1e-6*(10^d)*NA; % umol/dm^d -> mol/m^d -> spins/m^d

muB = 9.2740100783e-24; % Bohr magneton, J/T (CODATA 2018 value);
mu0 = 1.25663706212e-6; % magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
h = 6.62607015e-34; % reduced Planck constant, J/Hz (CODATA 2018)
ge = 2.00231930436256; % free-electron g factor (CODATA 2018 value)
hbar = h/2/pi;

D = (mu0/4/pi)*(muB*ge)^2/hbar; % m^3 s^-1

% Compute constants
if d==3
    c = -pi/2;
    Lam = 4/3/sqrt(3);
else
    c = cos(d*pi/6)*gamma(-d/3);
    fun = @(z) abs(1-3*z.^2).^(d/3);
    Lam = integral(fun,0,1);
end

% Compute background function
B = exp(4*pi/3*c*Lam*lambda*conc*D^(d/3)*abs(t*1e-6).^(d/3));
B = B(:);
output = B;

return

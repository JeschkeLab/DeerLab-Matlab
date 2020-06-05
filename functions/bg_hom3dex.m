%
% BG_HOM3DEX Excluded-volume model
%
%   info = BG_HOM3DEX
%   Returns an (info) table of model parameters and boundaries.
%
%   B = BG_HOM3DEX(t,param)
%   B = BG_HOM3DEX(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
%    ---------------------------------------------------------
%     Index  Parameter          Units   Lower   Upper   Start
%    ---------------------------------------------------------
%       1    Spin concentration  uM     0.01     1000    50  
%       2    Exclusion distance  nm     0.1       20     1 
%    ---------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = bg_hom3dex(t,param,lambda)

nParam = 2;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Spin concentration';
    info(1).Units = 'uM';
    info(1).Lower = 0.01;
    info(1).Upper = 1000;
    info(1).Start = 50;
    
    info(2).Index = 2;
    info(2).Parameter = 'Exclusion distance';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 20;
    info(2).Start = 1;
    
    output = struct2table(info);
    return
end

if nargin<3
    lambda = 1;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters (%d) does not match the number of model parameters (%d).',...
        length(param),nParam)
end

% Load precalculated reduction factor look-up table (Kattnig Eq.(18))
% To regenerate look-up table, use private/bg_hom3dex_alpha
persistent exvol
if isempty(exvol)
    load('bg_hom3dex','exvol');
end

% Get parameters
conc = param(1); % uM
R = param(2); % nm

NA = 6.02214076e23; % Avogadro constant, mol^-1
conc = conc*1e-6*1e3*NA; % umol/L -> mol/L -> mol/m^3 -> spins/m^3
ge = 2.00231930436256; % free-electron g factor (CODATA 2018 value)
mu0 = 1.25663706212e-6; % magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
muB = 9.2740100783e-24; % Bohr magneton, J/T (CODATA 2018 value);
h = 6.62607015e-34; % Planck constant, J/Hz (CODATA 2018)
hbar = h/2/pi;

A = (mu0/4/pi)*(ge*muB)^2/hbar; % Eq.(6); m^3 s^-1

% Calculate reduction factor (Eq.(18))
if R==0
    alpha = 1;
else
    dR = A*abs(t*1e-6)/(R*1e-9)^3; % unitless
    
    % Use interpolation of look-up table for small dR
    small = dR<max(exvol.dR);
    alpha = zeros(size(dR));
    alpha(small) = interp1(exvol.dR,exvol.alpha,dR(small),'makima');
    
    % For large dR, use limiting dR->inf expression
    alpha(~small) = 1 - (3/2/pi)*sqrt(3)./dR(~small);
end

K = 8*pi^2/9/sqrt(3)*A*abs(t*1e-6).*alpha; % Eq.(17)
B = exp(-lambda*conc*K); % Eq.(13)
B = B(:);

output = B;

return

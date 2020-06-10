function output = dd_randcoil(r,param)
%
% DD_RANDCOIL Random-coil model for an unfolded peptide/protein
%
%   info = DD_RANDCOIL
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_RANDCOIL(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the parameter array (param). The required parameters can also be found 
%   in the (info) structure. The end-to-end distance distribution is approximated 
%   by a Gaussian coil with proper mean distance, which is good for
%   sufficiently large N.
%
% PARAMETERS
%    -------------------------------------------------------------
%     Index  Parameter           Units   Lower    Upper    Start
%    -------------------------------------------------------------
%       1    Number of residues           2       1000      50
%       2    Segment length       nm      0.1     0.4       0.2
%       3    Scaling exponent             0.33    1         0.602
%    -------------------------------------------------------------
%
%   See: N. C. Fitzkee, G. D. Rose, PNAS 2004, 101(34), 12497-12502
%        https://doi.org/10.1073/pnas.0404236101
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

nParam = 3;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Number of residues';
    info(1).Units = '  ';
    info(1).Lower = 2;
    info(1).Upper = 1000;
    info(1).Start = 50;
    
    info(2).Index = 2;
    info(2).Parameter = 'Segment length';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 0.4;
    info(2).Start = 0.2;
    
    info(3).Index = 3;
    info(3).Parameter = 'Scaling exponent';
    info(3).Units = '  ';
    info(3).Lower = 0.33;
    info(3).Upper = 1;
    info(3).Start = 0.602;
    
    output = info;
    return
end
    
% Check that the number of parameters matches the model requirements
if numel(param)~=nParam
  error('The number of input parameters does not match the number of model parameters.');
end

% Parse input axis
if ~iscolumn(r)
  r = r.';
end
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

N = param(1); % number of residues
nu = param(2); % scaling exponent
R0 = param(3); % residue length

rsq = 6*(R0*N^nu)^2; % mean square end-to-end distance from radius of gyration
normFact = 3/(2*pi*rsq)^(3/2); % normalization prefactor
ShellSurf = 4*pi*r.^2; % spherical shell surface
Gaussian = exp(-3*r.^2/(2*rsq));
P = normFact*ShellSurf.*Gaussian;

% Normalize integral
if ~all(P==0)
P = P/trapz(r,P);
end
output = P;

end

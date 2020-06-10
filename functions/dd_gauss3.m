%
% DD_GAUSS3 Sum of three Gaussian distributions parametric model
%
%   info = DD_GAUSS3
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_GAUSS3(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    ------------------------------------------------------------------
%     Index  Parameter                Units   Lower    Upper    Start
%    ------------------------------------------------------------------
%       1    Center of 1st Gaussian     nm      1      20       2.5 
%       2    FWHM of 1st Gaussian       nm     0.2      5       0.5 
%       3    Amplitude of 1st Gaussian          0       1       0.3 
%       4    Center of 2nd Gaussian     nm      1      20       3.5 
%       5    FWHM of 2nd Gaussian       nm     0.2      5       0.5 
%       6    Amplitude of 2nd Gaussian          0       1       0.3 
%       7    Center of 3rd Gaussian     nm      1      20       5.0 
%       8    'FWHM of 3rd Gaussian      nm     0.2      5       0.5
%    ------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = dd_gauss3(r,param)

nParam = 8;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Center of 1st Gaussian';
    info(1).Units = 'nm';
    info(1).Lower = 1;
    info(1).Upper = 20;
    info(1).Start = 2.5;
    
    info(2).Index = 2;
    info(2).Parameter = 'FWHM of 1st Gaussian';
    info(2).Units = 'nm';
    info(2).Lower = 0.2;
    info(2).Upper = 5;
    info(2).Start = 0.5;
    
    info(3).Index = 3;
    info(3).Parameter = 'Amplitude of 1st Gaussian';
    info(3).Units = '  ';
    info(3).Lower = 0;
    info(3).Upper = 1;
    info(3).Start = 0.3;
    
    info(4).Index = 4;
    info(4).Parameter = 'Center of 2nd Gaussian';
    info(4).Units = 'nm';
    info(4).Lower = 1;
    info(4).Upper = 20;
    info(4).Start = 3.5;
    
    info(5).Index = 5;
    info(5).Parameter = 'FWHM of 2nd Gaussian';
    info(5).Units = 'nm';
    info(5).Lower = 0.2;
    info(5).Upper = 5;
    info(5).Start = 0.5;

    info(6).Index = 6;
    info(6).Parameter = 'Amplitude of 2nd Gaussian';
    info(6).Units = '  ';
    info(6).Lower = 0;
    info(6).Upper = 1;
    info(6).Start = 0.3;

    info(7).Index = 7;
    info(7).Parameter = 'Center of 3rd Gaussian';
    info(7).Units = 'nm';
    info(7).Lower = 1;
    info(7).Upper = 20;
    info(7).Start = 5.0;
    
    info(8).Index = 8;
    info(8).Parameter = 'FWHM of 3rd Gaussian';
    info(8).Units = 'nm';
    info(8).Lower = 0.2;
    info(8).Upper = 5;
    info(8).Start = 0.5;
    
    output = info;
    return
end
    
% Check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
fwhm = param([2 5 8]);
r0 = param([1 4 7]);
a = param([3 6]);
a(3) = max(1-sum(a),0);
P = multigaussfun(r,r0,fwhm,a);

output = P;

return

%
% DD_RICE4 Sum of four 3D-Rice distributions parametric model
%
%   info = DD_RICE4
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_RICE4(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% % PARAMETERS
%    -----------------------------------------------------------------
%     Index  Parameter                Units   Lower    Upper    Start
%    -----------------------------------------------------------------
%       1    Center of 1st Rician      nm       1       10      2.5 
%       2    Width of 1st Rician       nm      0.1       5      0.7 
%       3    Amplitude of 1st Rician            0        1      0.25 
%       4    Center of 2nd Rician      nm       1       10      3.5 
%       5    Width of 2nd Rician       nm      0.1       5      0.7 
%       6    Amplitude of 2nd Rician            0        1      0.25 
%       7    Center of 3rd Rician      nm       1       10      4.5 
%       8    Width of 3rd Rician       nm      0.1       5      0.7 
%       9    Amplitude of 3rd Rician            0        1      0.25 
%      10    Center of 4th Rician      nm       1       10      5.5 
%      11    Width of 4th Rician       nm      0.1     0.5      0.7 
%    -----------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = dd_rice4(r,param)

nParam = 11;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Center of 1st Rician';
    info(1).Units = 'nm';
    info(1).Lower = 1;
    info(1).Upper = 10;
    info(1).Start = 2.5;
    
    info(2).Index = 2;
    info(2).Parameter = 'Width of 1st Rician';
    info(2).Units = 'nm';
    info(2).Lower = 0.1;
    info(2).Upper = 5;
    info(2).Start = 0.7;
    
    info(3).Index = 3;
    info(3).Parameter = 'Amplitude of 1st Rician';
    info(3).Units = '  ';
    info(3).Lower = 0;
    info(3).Upper = 1;
    info(3).Start = 0.25;
    
    info(4).Index = 4;
    info(4).Parameter = 'Center of 2nd Rician';
    info(4).Units = 'nm';
    info(4).Lower = 1;
    info(4).Upper = 10;
    info(4).Start = 3.5;
    
    info(5).Index = 5;
    info(5).Parameter = 'Width of 2nd Rician';
    info(5).Units = 'nm';
    info(5).Lower = 0.1;
    info(5).Upper = 5;
    info(5).Start = 0.7;
    
    info(6).Index = 6;
    info(6).Parameter = 'Amplitude of 2nd Rician';
    info(6).Units = '  ';
    info(6).Lower = 0;
    info(6).Upper = 1;
    info(6).Start = 0.25;
    
    info(7).Index = 7;
    info(7).Parameter = 'Center of 3rd Rician';
    info(7).Units = 'nm';
    info(7).Lower = 1;
    info(7).Upper = 10;
    info(7).Start = 4.5;
    
    info(8).Index = 8;
    info(8).Parameter = 'Width of 3rd Rician';
    info(8).Units = 'nm';
    info(8).Lower = 0.1;
    info(8).Upper = 5;
    info(8).Start = 0.7;
    
    info(9).Index = 9;
    info(9).Parameter = 'Amplitude of 3rd Rician';
    info(9).Units = '  ';
    info(9).Lower = 0;
    info(9).Upper = 1;
    info(9).Start = 0.25;
    
    info(10).Index = 10;
    info(10).Parameter = 'Center of 4th Rician';
    info(10).Units = 'nm';
    info(10).Lower = 1;
    info(10).Upper = 10;
    info(10).Start = 5.5;
    
    info(11).Index = 11;
    info(11).Parameter = 'Width of 4th Rician';
    info(11).Units = 'nm';
    info(11).Lower = 0.1;
    info(11).Upper = 0.5;
    info(11).Start = 0.7;
    
    output = struct2table(info);
    return
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute non-central chi distribution with 3 degrees of freedom (a 3D Rician)
nu = param([1 4 7 10]);
sig = param([2 5 8 11]);
a = param([3 6 9]);
a(4) = max(1-sum(a),0);

P = multirice3d(r,nu,sig,a);

output = P;

return

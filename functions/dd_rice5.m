%
% DD_RICE5 Sum of five 3D-Rice distributions parametric model
%
%   info = DD_RICE5
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_RICE5(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
%    -----------------------------------------------------------------
%     Index  Parameter                Units   Lower    Upper    Start
%    -----------------------------------------------------------------
%       1    Center of 1st Rician       nm       1       10      2.5 
%       2    Width of 1st Rician        nm     0.1        5      0.7 
%       3    Amplitude of 1st Rician             0        1      0.2 
%       4    Center of 2nd Rician       nm       1       10      3.5 
%       5    Width of 2nd Rician        nm     0.1        5      0.7 
%       6    Amplitude of 2nd Rician             0        1      0.2 
%       7    Center of 3rd Rician       nm       1       10        4 
%       8    Width of 3rd Rician        nm     0.1        5      0.7 
%       9    Amplitude of 3rd Rician             0        1      0.2 
%      10    Center of 4th Rician       nm       1       10        5 
%      11    Width of 4th Rician        nm     0.1      0.5      0.7 
%      12    Amplitude of 4th Rician             0        1      0.2 
%      13    Center of 5th Rician       nm       1       10      5.5 
%      14    Width of 5th Rician        nm     0.1      0.5      0.7 
%    -----------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = dd_rice5(r,param)


nParam = 14;

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
    info(3).Start = 0.2;
    
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
    info(6).Start = 0.2;
    
    info(7).Index = 7;
    info(7).Parameter = 'Center of 3rd Rician';
    info(7).Units = 'nm';
    info(7).Lower = 1;
    info(7).Upper = 10;
    info(7).Start = 4;
    
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
    info(9).Start = 0.2;
    
    info(10).Index = 10;
    info(10).Parameter = 'Center of 4th Rician';
    info(10).Units = 'nm';
    info(10).Lower = 1;
    info(10).Upper = 10;
    info(10).Start = 5;
    
    info(11).Index = 11;
    info(11).Parameter = 'Width of 4th Rician';
    info(11).Units = 'nm';
    info(11).Lower = 0.1;
    info(11).Upper = 0.5;
    info(11).Start = 0.7;
        
    info(12).Index = 12;
    info(12).Parameter = 'Amplitude of 4th Rician';
    info(12).Units = '  ';
    info(12).Lower = 0;
    info(12).Upper = 1;
    info(12).Start = 0.2;
    
    info(13).Index = 13;
    info(13).Parameter = 'Center of 5th Rician';
    info(13).Units = 'nm';
    info(13).Lower = 1;
    info(13).Upper = 10;
    info(13).Start = 5.5;
    
    info(14).Index = 14;
    info(14).Parameter = 'Width of 5th Rician';
    info(14).Units = 'nm';
    info(14).Lower = 0.1;
    info(14).Upper = 0.5;
    info(14).Start = 0.7;
    
    output = struct2table(info);
    return
    
    % If no inputs given, return info about the parametric model
    info.model  = 'Five 3D-Rice distributions';
    info.nparam  = nParam;
    info.parameters(1).name = ['Center ',char(957),'1 1st component'];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Spread ',char(963),'1 1st component'];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
        
    info.parameters(3).name = 'Relative amplitude a1 1st component';
    info.parameters(3).range = [0 1];
    info.parameters(3).default = 0.2;
    info.parameters(3).units = '';

    info.parameters(4).name = ['Center ',char(957),'2 2nd component'];
    info.parameters(4).range = [1 10];
    info.parameters(4).default = 3.5;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = ['Spread ',char(963),'2 2nd component'];
    info.parameters(5).range = [0.1 5];
    info.parameters(5).default = 0.7;
    info.parameters(5).units = 'nm';
        
    info.parameters(6).name = 'Relative amplitude a2 2nd component';
    info.parameters(6).range = [0 1];
    info.parameters(6).default = 0.2;
    info.parameters(6).units = '';
    
    info.parameters(7).name = ['Center ',char(957),'3 3rd component'];
    info.parameters(7).range = [1 10];
    info.parameters(7).default = 4.0;
    info.parameters(7).units = 'nm';
    
    info.parameters(8).name = ['Spread ',char(963),'3 3rd component'];
    info.parameters(8).range = [0.1 5];
    info.parameters(8).default = 0.7;
    info.parameters(8).units = 'nm';
         
    info.parameters(9).name = 'Relative amplitude a3 3rd component';
    info.parameters(9).range = [0 1];
    info.parameters(9).default = 0.2;
    info.parameters(9).units = '';
    
    info.parameters(10).name = ['Center ',char(957),'4 4th component'];
    info.parameters(10).range = [1 10];
    info.parameters(10).default = 5;
    info.parameters(10).units = 'nm';
    
    info.parameters(11).name = ['Spread ',char(963),'4 4th component'];
    info.parameters(11).range = [0.1 5];
    info.parameters(11).default = 0.7;
    info.parameters(11).units = 'nm';
     
    info.parameters(12).name = 'Relative amplitude a4 4th component';
    info.parameters(12).range = [0 1];
    info.parameters(12).default = 0.2;
    info.parameters(12).units = '';
    
    info.parameters(13).name = ['Center ',char(957),'5 5th component'];
    info.parameters(13).range = [1 10];
    info.parameters(13).default = 5.5;
    info.parameters(13).units = 'nm';
    
    info.parameters(14).name = ['Spread ',char(963),'5 5th component'];
    info.parameters(14).range = [0.1 5];
    info.parameters(14).default = 0.7;
    info.parameters(14).units = 'nm';
        
    output = info;
    return;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute non-central chi distribution with 3 degrees of freedom (a 3D Rician)
nu = param([1 4 7 10 13]);
sig = param([2 5 8 11 14]);
a = param([3 6 9 12]);
a(5) = max(1-sum(a),0);

P = multirice3d(r,nu,sig,a);

output = P;

return

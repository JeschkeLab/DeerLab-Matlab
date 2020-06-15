%
% DD_GAUSS5 Sum of five Gaussian distributions parametric model
%
%   info = DD_GAUSS5
%   Returns an (info) table of model parameters and boundaries.
%
%   P = DD_GAUSS5(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
%    ------------------------------------------------------------------
%     Index  Parameter                  Units   Lower   Upper    Start
%    ------------------------------------------------------------------
%       1    Center of 1st Gaussian      nm      1       20       2.5 
%       2    FWHM of 1st Gaussian        nm     0.2       5       0.5 
%       3    Amplitude of 1st Gaussian           0        1       0.2 
%       4    Center of 2nd Gaussian      nm      1       20        3 
%       5    FWHM of 2nd Gaussian        nm     0.2       5       0.5 
%       6    Amplitude of 2nd Gaussian           0        1       0.2 
%       7    Center of 3rd Gaussian      nm      1       20       3.5 
%       8    FWHM of 3rd Gaussian        nm     0.2       5       0.5 
%       9    Amplitude of 3rd Gaussian           0        1       0.2 
%      10    Center of 4th Gaussian      nm      1       20        4 
%      11    FWHM of 4th Gaussian        nm     0.2       5       0.5 
%      12    Amplitude of 4th Gaussian           0        1       0.2 
%      13    Center of 5th Gaussian      nm      1       20        5 
%      14    FWHM of 5th Gaussian        nm     0.2       5       0.5 
%      15    Amplitude of 5th Gaussian           0        1       0.2 
%    ------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_gauss5(r,param)

nParam = 15;

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
    info(3).Start = 0.2;
    
    info(4).Index = 4;
    info(4).Parameter = 'Center of 2nd Gaussian';
    info(4).Units = 'nm';
    info(4).Lower = 1;
    info(4).Upper = 20;
    info(4).Start = 3.0;
    
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
    info(6).Start = 0.2;

    info(7).Index = 7;
    info(7).Parameter = 'Center of 3rd Gaussian';
    info(7).Units = 'nm';
    info(7).Lower = 1;
    info(7).Upper = 20;
    info(7).Start = 3.5;
    
    info(8).Index = 8;
    info(8).Parameter = 'FWHM of 3rd Gaussian';
    info(8).Units = 'nm';
    info(8).Lower = 0.2;
    info(8).Upper = 5;
    info(8).Start = 0.5;
    
    info(9).Index = 9;
    info(9).Parameter = 'Amplitude of 3rd Gaussian';
    info(9).Units = '  ';
    info(9).Lower = 0;
    info(9).Upper = 1;
    info(9).Start = 0.2;    
    
    info(10).Index = 10;
    info(10).Parameter = 'Center of 4th Gaussian';
    info(10).Units = 'nm';
    info(10).Lower = 1;
    info(10).Upper = 20;
    info(10).Start = 4.0;
    
    info(11).Index = 11;
    info(11).Parameter = 'FWHM of 4th Gaussian';
    info(11).Units = 'nm';
    info(11).Lower = 0.2;
    info(11).Upper = 5;
    info(11).Start = 0.5;

    info(12).Index = 12;
    info(12).Parameter = 'Amplitude of 4th Gaussian';
    info(12).Units = '  ';
    info(12).Lower = 0;
    info(12).Upper = 1;
    info(12).Start = 0.2; 
    
    info(13).Index = 13;
    info(13).Parameter = 'Center of 5th Gaussian';
    info(13).Units = 'nm';
    info(13).Lower = 1;
    info(13).Upper = 20;
    info(13).Start = 5.0;
    
    info(14).Index = 14;
    info(14).Parameter = 'FWHM of 5th Gaussian';
    info(14).Units = 'nm';
    info(14).Lower = 0.2;
    info(14).Upper = 5;
    info(14).Start = 0.5;

    info(15).Index = 15;
    info(15).Parameter = 'Amplitude of 5th Gaussian';
    info(15).Units = '  ';
    info(15).Lower = 0;
    info(15).Upper = 1;
    info(15).Start = 0.2; 
    
    output = info;
    return

end

if nargin~=2
    error('Model requires two input arguments.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
fwhm = param([2 5 8 11 14]);
r0 = param([1 4 7 10 13]);
a = param([3 6 9 12 15]);
P = multigaussfun(r,r0,fwhm,a);

output = P;

return

%
% EX_5PDEER 7-pulse DEER experiment model 
%
%   info = EX_7PDEER()
%   Returns an (info) table of model parameters and boundaries.
%
%   pathways = EX_7PDEER(param)
%   Computes the dipolar pathway information array according to the paramater
%   array (param).
%
%
% PARAMETERS
%    ----------------------------------------------------------------------------
%     Index  Parameter                               Units   Lower  Upper  Start
%    ----------------------------------------------------------------------------
%       1    Amplitude of unmodulated components               0      1     0.3 
%       2    Amplitude of 1st modulated pathway                0      1     0.5 
%       3    Amplitude of 2nd modulated pathway                0      1     0.3 
%       4    Amplitude of 3rd modulated pathway                0      1     0.2 
%       5    Refocusing time of 2nd modulated pathway  us      0      20    1.5 
%       6    Refocusing time of 3rd modulated pathway  us      0      20    3.5 
%    ----------------------------------------------------------------------------
%       

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = ex_7pdeer(param)

nParam = 6;

if nargin>1
    error('Model requires one input argument.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Amplitude of unmodulated components';
    info(1).Units = '  ';
    info(1).Lower = 0;
    info(1).Upper = 1;
    info(1).Start = 0.3;
    
    info(2).Index = 2;
    info(2).Parameter = 'Amplitude of 1st modulated pathway';
    info(2).Units = '  ';
    info(2).Lower = 0;
    info(2).Upper = 1;
    info(2).Start = 0.5;
    
    info(3).Index = 3;
    info(3).Parameter = 'Amplitude of 2nd modulated pathway';
    info(3).Units = '  ';
    info(3).Lower = 0;
    info(3).Upper = 1;
    info(3).Start = 0.3;

    info(4).Index = 4;
    info(4).Parameter = 'Amplitude of 3rd modulated pathway';
    info(4).Units = '  ';
    info(4).Lower = 0;
    info(4).Upper = 1;
    info(4).Start = 0.2;
    
    info(5).Index = 5;
    info(5).Parameter = 'Refocusing time of 2nd modulated pathway';
    info(5).Units = 'us';
    info(5).Lower = 0;
    info(5).Upper = 20;
    info(5).Start = 1.5;

    info(6).Index = 6;
    info(6).Parameter = 'Refocusing time of 3rd modulated pathway';
    info(6).Units = 'us';
    info(6).Lower = 0;
    info(6).Upper = 20;
    info(6).Start = 3.5;
    
    output = struct2table(info);
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end
param = param(:).';

% Extract parameters
lambda = param(1:4);
lambda = lambda/sum(lambda);
T0 = [0 param(5:6)];

% Dipolar pathways
pathways(1,:) = [lambda(1) NaN];
pathways(2,:) = [lambda(2) T0(1)];
pathways(3,:) = [lambda(3) T0(2)];
pathways(4,:) = [lambda(4) T0(3)];
output = pathways;

end

%
% EX_5PDEER 5-pulse DEER experiment model
%
%   info = EX_5PDEER(t)
%   Returns an (info) structure containing the specifics of the model, including
%   a list of parameters.
%
%   pathways = EX_5PDEER(t,param)
%   Computes the dipolar pathway information array according to the paramater
%   array (param).
%
%
% PARAMETERS
%    ----------------------------------------------------------------------------------------
%     Index  Parameter                                 Units    Lower    Upper     Start
%    ----------------------------------------------------------------------------------------
%       1    Amplitude of unmodulated components                   0        1        0.3
%       2    Amplitude of 1st modulated pathway                    0        1        0.3
%       3    Amplitude of 2nd modulated pathway                    0        1        0.3
%       4    Refocusing time of 2nd modulated pathway   us  max(t)/2-2  max(t)/2+2  max(t)/2
%    ----------------------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = ex_5pdeer(t,param)

nParam = 4;

if nargin>2
    error('Model requires one or two input arguments.')
end

if nargin==1
    % If no inputs given, return info about the parametric model
    info(1).Index = 1;
    info(1).Parameter = 'Amplitude of unmodulated components';
    info(1).Units = '  ';
    info(1).Lower = 0;
    info(1).Upper = 1;
    info(1).Start = 0.4;
    
    info(2).Index = 2;
    info(2).Parameter = 'Amplitude of 1st modulated pathway';
    info(2).Units = '  ';
    info(2).Lower = 0;
    info(2).Upper = 1;
    info(2).Start = 0.4;
    
    info(3).Index = 3;
    info(3).Parameter = 'Amplitude of 2nd modulated pathway';
    info(3).Units = '  ';
    info(3).Lower = 0;
    info(3).Upper = 1;
    info(3).Start = 0.2;
    
    info(4).Index = 4;
    info(4).Parameter = 'Refocusing time of 2nd modulated pathway';
    info(4).Units = 'us';
    info(4).Lower = max(t)/2 - 2;
    info(4).Upper = max(t)/2 + 2;
    info(4).Start = max(t)/2;

    output = struct2table(info);
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Extract parameter
lambda = param(1:3);
lambda = lambda/sum(lambda);
T0 = [0 param(4)];

% Dipolar pathways
pathways(1,:) = [lambda(1) NaN];
pathways(2,:) = [lambda(2) T0(1)];
pathways(3,:) = [lambda(3) T0(2)];
output = pathways;

end

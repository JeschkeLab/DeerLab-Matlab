%
% APT Approximate Pake Transformation of dipolar spectroscopy data
%
%   P = APT(S,K)
%   Computes the distance distribution corresponding to the signal (S)
%   according to the APT kernel (K). The kernel must be a structure as
%   generated by the aptkernel() function.
%
%   P = APT(S,K,DDS)
%   Specify the distance-domain smoothing parameter (DDS) by passing as a 
%   scalar value.
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [Distribution,UniformDistanceAxis] = apt(Signal,APTkernel,DistDomainSmoothing)

%--------------------------------------------------------------------------
% Parse & Validate Required Input
%--------------------------------------------------------------------------
if ~isa(APTkernel,'struct')
    error('The input APTkernel must be a a valid structure.')
end

if iscolumn(Signal)
   Signal = Signal'; 
end
validateattributes(Signal,{'numeric'},{'2d','nonempty'})

if nargin<2 || isempty(DistDomainSmoothing)
    DistDomainSmoothing = 0.05;
else
    validateattributes(DistDomainSmoothing,{'numeric'},{'scalar','nonnegative'})
end

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({Signal,APTkernel,DistDomainSmoothing});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [Distribution,UniformDistanceAxis] = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
% APT algorithm
%--------------------------------------------------------------------------

%Get APT kernel data
Kernel = APTkernel.Base;
NormConstant = APTkernel.NormalizationFactor;
APT_FrequencyAxis = APTkernel.FreqAxis;
APT_TimeAxis = APTkernel.TimeAxis;
Crosstalk = APTkernel.Crosstalk;

%Compute frequency distribution
[FreqDimension,~] = size(Kernel);
FreqDistribution = zeros(1,FreqDimension);
for k=1:FreqDimension
    FreqDistribution(k) = FreqDistribution(k)+sum(Kernel(k,:).*Signal.*APT_TimeAxis)/NormConstant(k);
end

%Perform crosstalk correction
APTdistribution = Crosstalk\FreqDistribution'; % crosstalk correction, eqn [22]

%Map fequencies to distances
Freq2Dist = zeros(length(APT_FrequencyAxis),1); % initialize distance axis (mapping of dipolar frequencies to distances)
for k = 1:length(Freq2Dist)
    Freq2Dist(k) = (52.04/APT_FrequencyAxis(k))^(1/3);
end
MappedDistances = APT_FrequencyAxis/2; % initialize distance axis (mapping of dipolar frequencies to distances)
for k = 1:length(MappedDistances)
    MappedDistances(k) = (52.04/MappedDistances(k))^(1/3);
end
APTdistribution = interp1(Freq2Dist,APTdistribution,MappedDistances,'pchip',0);
APTdistribution = APTdistribution';

%Perform distance-domain smoothing filtering
FilteredAPTdistribution = zeros(length(APTdistribution),1);
for k = 1:length(APTdistribution)
    DDSfilter = (MappedDistances - MappedDistances(k)*ones(1,length(MappedDistances)))/DistDomainSmoothing;
    DDSfilter = exp(-DDSfilter.*DDSfilter);
    DDSfilter = DDSfilter';
    FilteredAPTdistribution(k) = sum(DDSfilter.*APTdistribution)/sum(DDSfilter);
end

%Normalize integral
for k=1:length(APTdistribution)
    FilteredAPTdistribution(k) = FilteredAPTdistribution(k)/(MappedDistances(k))^4;
end
%Renormalize
FilteredAPTdistribution = FilteredAPTdistribution/max(FilteredAPTdistribution);

%Interpolate to uniform distance axis
UniformDistanceAxis = linspace(min(MappedDistances),max(MappedDistances),length(Signal));
Distribution = uniformgrain(MappedDistances,FilteredAPTdistribution,UniformDistanceAxis);

%Normalize to unity integral
Distribution = Distribution/sum(Distribution)/mean(diff(UniformDistanceAxis));

%Make the distribution a column
Distribution = Distribution';

%Store output result in the cache
Output = {Distribution,UniformDistanceAxis};
cachedData = addcache(cachedData,hashKey,Output);


end


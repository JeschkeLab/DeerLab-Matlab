%
% CORRECTZEROTIME Zero-time correction of dipolar spectroscopy signals
%
%   tc = CORRECTZEROTIME(V,t)
%   Determines the zero time of a dipolar signal (V) and corrects
%   the time axis (t) for it, returning a corrected time axis (tc).
%   If t is in ns/us, tc will be be in ns/us as well.
%
%   [tc,t0,pos] = CORRECTZEROTIME(V,t)
%   Returns the corrected time axis (tc), the zero-time (t0) in ns/us and the 
%   zero-time array index (pos). 
%
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [tcorr,t0,idxt0] = correctzerotime(V,t)

if nargin<2
    error('Not enough input arguments.')
end
if any(t<0)
   t = t + abs(min(t)); 
end
validateattributes(V,{'numeric'},{'2d'},mfilename,'S')
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t')

%Use column vectors
t = t(:);
V = V(:);

%Generate finely-grained interpolated signal and time axis
resolution = 4;
tfine = linspace(min(t),max(t),(numel(t)-1)*resolution+1);
Vfine = interp1(t,real(V),tfine,'spline');

% Determine maximum
[~,maxPos] = max(Vfine);
idxt0 = 1;
%If maximum is not the first or last point, then do moment analysis
if maxPos>1 && maxPos<length(Vfine)
    % Determine the width of the interval to use for the integral
    if maxPos<length(Vfine)/2
        maxDelta = floor((maxPos-1)/2);
    else
        maxDelta = floor((length(Vfine)-maxPos)/2);
    end
    offsetrange = -maxDelta:maxDelta;
    
    %Look around the maximum for lowest moment integral
    minIntegral = realmax;
    idxt0 = maxPos;
    for idx = maxPos-maxDelta:maxPos+maxDelta
        %Query new candidate for zero position and integrate
        dt = tfine(idx+offsetrange)-tfine(idx);
        Integral = sum(Vfine(idx+offsetrange).*dt);
        %If integral is lower than prior best, then update new candidate
        if abs(Integral) < minIntegral
            minIntegral = abs(Integral);
            idxt0 = idx;
        end
    end
end
t0 = tfine(idxt0);


% Correct time axis
tcorr = t - t0;
[~,idxt0] = min(abs(tcorr));

end

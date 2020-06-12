%
% CORRECTZEROTIME Zero-time correction of dipolar spectroscopy signals
%
%   tc = CORRECTZEROTIME(V,t)
%   Determines the zero time of a dipolar signal (V) and corrects
%   the time axis (t) for it, returning a corrected time axis (tc).
%
%   [tc,t0] = CORRECTZEROTIME(V,t)
%   [tc,t0,pos] = CORRECTZEROTIME(V,t)
%   Returns the corrected time axis (tc), the zero-time (t0) in ns/us and the
%   array index (pos) that is closest to the zero time.
%
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [tcorr,t0,idxt0] = correctzerotime(V,t)

% Validate input arguments
%-------------------------------------------------------------------------------
if nargin<2
    error('At least two inputs are required: V (signal) and t (time axis).')
elseif nargin>2
    error('At most two inputs are possible: V (signal) and t (time axis).');
end
validateattributes(V,{'numeric'},{'vector'},mfilename,'V')
validateattributes(t,{'numeric'},{'vector'},mfilename,'t')


% Determine zero-time via first-moment integral
%-------------------------------------------------------------------------------
% Use column vectors
V = V(:);
t = t(:);

% Generate finely-grained interpolated signal and time axis
resolution = 10;
t_ = linspace(min(t),max(t),(numel(t)-1)*resolution+1);
V_ = interp1(t,real(V),t_,'spline');

% Determine location of maximum
[~,idxmax] = max(V_);

if idxmax~=1 && idxmax~=numel(V)
    % If maximum is not the first or last point, then do moment analysis
    % (i.e. minimize the magntitude of the first-moment integral)
    
    % Determine the interval to use for the integral
    maxDelta = floor(min(idxmax-1,length(V_)-idxmax)/2);
    offsetrange = -maxDelta:maxDelta;
    
    % Look around the maximum for lowest-magnitude moment integral
    minIntegral = inf;
    absInt = @(i)abs(sum(V_(offsetrange+i).*(t_(offsetrange+i)-t_(i))));
    % Run over all zero-position candidates and evalulate the integral
    for idx = idxmax-maxDelta:idxmax+maxDelta
        absIntegral = absInt(idx);
        % If |integral| is smaller than prior best, then update candidate
        if absIntegral < minIntegral
            minIntegral = absIntegral;
            idxt0 = idx;
        end
    end
    
else
    % Return maximum if it's the first or the last point
    idxt0 = idxmax;
    
end

t0 = t_(idxt0);

% Correct time axis
tcorr = t - t0;

% Determine index of time axis point closest to t0
[~,idxt0] = min(abs(tcorr));

end

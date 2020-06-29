%
% PARSEOPTIONAL Parser of property-value pair inputs
%
%   varargout = PARSEOPTIONAL(props,varargin)
%   Takes all inputs in (varargin) and identifies the allowed properties
%   specified by the cell array (props).
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function varargout = parseoptions(Names,varargin)

% If only one property passed put it into a cell array
if length(Names)==1 && ~iscell(Names)
    Names = {Names};
end

% Ensure that all properties are lower case
Names = lower(Names);

% Pre-allocate output container
varargout = cell(length(Names),1);

% If no varargins passed, just return
if isempty(varargin{1}), return; end

if length(varargin)==1
    varargin = varargin{1};
    if all(cellfun(@length,varargin)>1) && all(cellfun(@(x)isa(x,'cell'),varargin))
        varargin = varargin{1};
    end
end
if isempty(varargin{1})
    varargout = cell(nargout,1);
    return
end

% Parse property-value pairs in varargin
% ------------------------------------------------------------------------------
for i = 1:2:length(varargin(:))
    currentName = varargin{i};
    if isa(currentName,'char')
        % If property does exist then add the value to output list
        if any(strcmpi(currentName,Names))
            argoutidx = strcmpi(currentName,Names);
            varargout{argoutidx} = varargin{i+1};
        end
    else
        % If the varargin does not contain Name-alue pairs, then show error
        error('DeerLab:parseoptions','The input is not a valid name-value pair.')
    end
    
end

end

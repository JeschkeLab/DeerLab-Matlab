
function varargout = parseoptions(options,calls)

if nargin<2
    calls = [];
end

% Get environment information
DeerLab_path = which('dipolarkernel.m');
functions_path = DeerLab_path(1:end-numel('dipolarkernel.m'));

%Get hardcoded list of options accepted by DeerLab functions
optionslist = DeerLabOptions();

% Get names of all DeerLab functions with options
DeerLabfcns = fieldnames(optionslist);

% Check the stack of calls
stack = dbstack;
% Get name of DeerLab function currently calling checkoptions()...
caller = stack(2).name;
% Remove subfunction name if parseoptions is called from a subfunction
if contains(caller,'/')
    idx = strfind(caller,'/');
    caller = caller(1:idx-1);
end
% ... as well as other functions that have been called
called = [{stack.name}];
calledfcns = [];
calledfcns = [calledfcns calls];
for i=1:numel(called)
    % See if called function is a DeerLab function
    idx = find(strcmpi(DeerLabfcns,called{i}));
    % Store the name of the called DeerLab function
    calledfcns = [calledfcns DeerLabfcns(idx)];
    % and remove from the list of available functions
    DeerLabfcns(idx) = [];
end


cache = load(fullfile(functions_path,'private','dependencies_cache.mat'));
nestedfcns = cache.dependencies.(caller);
calledfcns = [calledfcns nestedfcns];

% Check that the caller function is a DeerLab function
if all(strcmpi(DeerLabfcns,caller))
    error('The caller function ''%s'' is not a valid DeerLab function.',caller)
end

% Convert name-value pairs into struct
if iscell(options)
    Properties = optionslist.(caller);
    options = pairs2struct(options);
end

if ~isempty(options)
fields = fieldnames(options);
found = false;
for i=1:numel(fields)
    
    %Check if the options was used in a called function
    for j=1:numel(calledfcns)
        allowedoptions = optionslist.(calledfcns{j});
        found(j) =  any(strcmpi(fields{i},allowedoptions));
    end
    found = any(found);
    
    
    if ~found
        % If not search in all other functions
        foundfcn = [];
        found = false;
        for j=1:numel(DeerLabfcns)
            allowedoptions = optionslist.(DeerLabfcns{j});
            %If function allows the option then note that the option was found
            if any(strcmpi(fields{i},allowedoptions))
                if ~isempty(foundfcn)
                    foundfcn = [foundfcn ',''' DeerLabfcns{j} ''''];
                else
                    foundfcn = ['''' DeerLabfcns{j} ''''];
                end
                found = true;
            end
        end
        
        if found
            %If the option was found in some other function then recommend the
            %correct function
            error('''%s'' is not an option of the function ''%s'', but of the function(s) %s.',fields{i},caller,foundfcn)
        else
            %Otherwise prompt error that option is not valid
            error('''%s'' is not a valid DeerLab function option.',fields{i});
        end
    end
    
end
end

% Parse the options and return outputs
outputoptions = optionslist.(caller);
varargout = cell(numel(outputoptions),1);
for i=1:numel(outputoptions)
    
    if isfield(options,lower(outputoptions{i}))
        varargout{i} = options.(lower(outputoptions{i}));
    else
        varargout{i} = [];
    end
end

% Name-value pairs -> strucure convertor
%-------------------------------------------
    function optionsstruct = pairs2struct(optionspairs)
        
        if length(optionspairs)==1 && iscell(optionspairs)
            optionspairs = optionspairs{1};
        end
        
        optionsstruct = [];
        for ii = 1:2:length(optionspairs(:))
            currentProperty = lower(optionspairs{ii});
            if isa(currentProperty,'char')
                optionsstruct.(currentProperty) = optionspairs{ii+1};
            end
        end
    end


% List of DeerLab functions options
%-------------------------------------------
    function optionslist = DeerLabOptions()
        
        %==================================================
        %                   fitparamodel()
        %==================================================
        optionslist.fitparamodel = {'Solver','Algorithm','MaxIter','Verbose','MaxFunEvals',...
                                    'TolFun','GlobalWeights','MultiStart','Rescale'};
        
        %==================================================
        %                   fitregmodel()
        %==================================================
        optionslist.fitregmodel = {'TolFun','Solver','NonNegConstrained','Verbose','MaxFunEvals',...
                                   'MaxIter','HuberParam','GlobalWeights','OBIR','RegOrder',...
                                   'NormP'};
        
        %==================================================
        %                   selregparam()
        %==================================================
        optionslist.selregparam = {'TolFun','NonNegConstrained','NoiseLevel','GlobalWeights',...
                                   'HuberParameter','Range','RegOrder','Search'};
        
        %==================================================
        %                   dipolarkernel()
        %==================================================
        optionslist.dipolarkernel = {'ExcitationBandwidth','OvertoneCoeffs','g','Method','nKnots',...
                                     'Cache','Renormalize'};
                                 
        %==================================================
        %                   aptkernel()
        %==================================================
        optionslist.aptkernel = {'ExcitationBandwidth'};
                   
        %==================================================
        %                   backgroundstart()
        %==================================================
        
        optionslist.backgroundstart = {'SearchStart','SearchEnd','EndCutOff'};
                   
        %==================================================
        %                   dipolarbackground()
        %==================================================
        
        optionslist.dipolarbackground = {'OvertoneCoeffs','Renormalize'};
                   
        %==================================================
        %                   dipolarsignal()
        %==================================================
        
        optionslist.dipolarsignal = {'NoiseLevel','g','Scale','Overtones','Phase'};
                   
        %==================================================
        %                   fitbackground()
        %==================================================
        
        optionslist.fitbackground = {'LogFit','InitialGuess','ModDepth','Solver'};
                   
        %==================================================
        %                   fitmultimodel()
        %==================================================
        
        optionslist.fitmultimodel = {'Upper','Lower','Background','internal::parselater'};
                   
        %==================================================
        %                   fitsignal()
        %==================================================
        
        optionslist.fitsignal = {'RegParam','RegType','alphaOptThreshold','TolFun','normP','GlobalWeights'};
                   
        %==================================================
        %                   longpass()
        %==================================================
        
        optionslist.longpass = {};
                   
        %==================================================
        %                   selectmodel()
        %==================================================
        
        optionslist.selectmodel = {'GlobalWeights'};
                   
        %==================================================
        %                   winlowpass()
        %==================================================
        
        optionslist.winlowpass = {'MinimalAttenuation','ForwardBackward'};
                   
        %==================================================
        %                   bootan()
        %==================================================
        
        optionslist.bootan = {'Verbose','Resampling'};
                   
        %==================================================
        %                   fftspec()
        %==================================================
        
        optionslist.fftspec = {'Type','ZeroFilling','Apodization'};
                   
        %==================================================
        %                   regparamrange()
        %==================================================
        
        optionslist.regparamrange = {'NoiseLevel','Resolution'};
                   
        %==================================================
        %                   prepvalidation()
        %==================================================
        
        optionslist.prepvalidation = {'randperm'};
                   
        %==================================================
        %                   sensitivan()
        %==================================================
        
        optionslist.sensitivan = {'AxisHandle','RandPerm','dynamicStats','Verbose'};
        
        %==================================================
        %                   snlls()
        %==================================================
        
        optionslist.snlls = {'alphaOptThreshold','RegOrder','RegParam','RegType',...
                             'nonLinSolver','LinSolver','forcepenalty',...
                             'nonLinMaxIter','nonLinTolFun','LinMaxIter','LinTolFun',...
                             'multiStarts'};
    end




% Cache regenerator of DeerLab function dependencies
%----------------------------------------------------
    function regenDependencies()
        %
        % This function regenerates the file 'dependencies_cache.mat' which
        % contains a structure with all the DeerLab function dependencies required
        % for the option system to work.
        %
        % This cache contains the information from the MATLAB API function
        % matlab.codetools.requiredFilesAndProducts which is very slow (19.5s) in
        % comparison to the cached version (0.01s).
        %
        % Must be called by runnning the function manually.
        
        % Get list of DeerLab functions with options
        optionslist = DeerLabOptions();
        files = fieldnames(optionslist);
        % Use MATLAB API to get the function dependencies

        for ii=1:numel(files)
            % Get dependencies
            list = matlab.codetools.requiredFilesAndProducts(files{ii});
            % Prune absolute paths and keep just function names
            [~,list] = cellfun(@(c)fileparts(c),list,'UniformOutput',false);

            fcnlist = [];
            for jj=1:numel(list)
               if any(strcmp(list{jj},files))
                  fcnlist = [fcnlist list(jj)]; 
               end
            end
            % Store in structure
            dependencies.(files{ii}) = fcnlist;
        end
        
        % Store the dependencies into a file
        save(fullfile(functions_path,'private','dependencies_cache.mat'),'dependencies')
        
    end



end

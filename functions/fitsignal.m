%
% FITSIGNAL  Fit model to dipolar time-domain trace
%
%   [Vfit,Pfit,Bfit,parfit,modfitci,parci,stats] = FITSIGNAL(V,t,r,dd,bg,ex,par0,lb,ub)
%   __ = FITSIGNAL(V,t,r,dd,bg,ex,par0)
%   __ = FITSIGNAL(V,t,r,dd,bg,ex)
%   __ = FITSIGNAL(V,t,r,dd,bg)
%   __ = FITSIGNAL(V,t,r,dd)
%   __ = FITSIGNAL(V,t,r)
%   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__},par0,lb,ub)
%   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__},par0)
%   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},{ex1,ex2,__})
%   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd,{bg1,bg2,__},ex)
%   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r,dd)
%   __ = FITSIGNAL({V1,V2,__},{t1,t2,__},r)
%   __ = FITSIGNAL(___,'Property',Values,___)
%
%   Fits a dipolar model to the experimental signal V with time axis t, using
%   distance axis r. The model is specified by the distance distribution (dd),
%   the background (bg), and the experiment (ex).
%
%   If multiple signals (V1,V2,...) and their corresponding time axes (t1,t2,...)
%   are given, they will be fitted globally with a single distance distribution (dd).
%   For each signal, a specific background (bg1,bg2,...) and experiment (ex1,ex2)
%   models can be assigned.
%
%   FITSIGNAL can handle both parametric and non-parametric distance
%   distribution models.
%
%  Input:
%    V      time-domain signal to fit (N-element vector)
%    t      time axis, in microseconds (N-element vector)
%    r      distance axis, in nanometers (M-element vector)
%    dd     distance distribution model (default 'P')
%           - 'P' to indicate non-parametric distribution
%           - function handle to parametric distribution model (e.g. @dd_gauss)
%           - 'none' to indicate no distribution, i.e. only background
%    bg     background model (default @bg_hom3d)
%           - function handle to parametric background model (e.g. @bg_hom3d)
%           - 'none' to indicate no background decay
%    ex     experiment model (default @ex_4pdeer)
%           - function handle to experiment model (e.g. @ex_4pdeer)
%           - 'none' to indicate simple dipolar oscillation (mod.depth = 1)
%    par0   starting parameters, 3-element cell array {par0_dd,par0_bg,par0_ex}
%           default: {[],[],[]} (automatic choice)
%    lb     lower bounds for parameters, 3-element cell array (lb_dd,lb_bg,lb_ex)
%           default: {[],[],[]} (automatic choice)
%    ub     upper bounds for parameters, 3-element cell array (ub_dd,ub_bg,ub_ex)
%           default: {[],[],[]} (automatic choice)
%
%  Output:
%    Vfit     fitted time-domain signal
%    Pfit     fitted distance distribution
%    Bfit     fitted background decay
%    parfit   structure with fitted parameters
%                .dd  fitted parameters for distance distribution model
%                .bg  fitted parameters for background model
%                .ex  fitted parameters for experiment model
%    modfitci structure with confidence intervals for Vfit, Bfit and Pfit
%    parci    structure with confidence intervals for parameter, similar to parfit
%    stats    goodness of fit statistical estimators, N-element structure array
%
%  Name-value pairs:
%
%   'Rescale'       - Enable/Disable optimization of the signal scale
%   'normP'         - Enable/Disable re-normalization of the fitted distribution
%   'Display'       - Enable/Disable plotting and printing of results
%   See "help snnls" for more options.
%
% Example:
%    Vfit = fitsignal(Vexp,t,r,@dd_gauss,@bg_hom3d,@ex_4pdeer)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [Vfit,Pfit,Bfit,parfit,modfituq,paruq,stats] = fitsignal(Vexp,t,r,varargin)

if nargin<3
    error('At least three inputs (V,t,r) must be specified.');
end

% Parse inputs in the varargin
[dd_model,bg_model,ex_model,par0,lb,ub,options] = parseinputs(varargin);

% Validate the V,t,r inputs and check global fit settings
nSignals = validateVtr;

% Optional arguments
% ====================

% Default optional settings
regtype = 'tikh';
regparam = 'aic';
DisplayResults = nargout==0;
normP = true;

% Parse and validate options passed by the user, if the user has specified
% any, this call will overwrite the defaults above
parsevalidate(options)

calculateCI = nargout>=5 || DisplayResults;
computeStats = nargout>4 || DisplayResults;
verbose = 'off';

% Get information about distance distribution parameters
% ======================================================
parametricDistribution  = isa(dd_model,'function_handle');
parfreeDistribution = strcmp(dd_model,'P');
includeForeground = ~strcmp(dd_model,'none');

[par0_dd,lower_dd,upper_dd,N_dd] = getmodelparams(dd_model);

if includeForeground && ~parametricDistribution && ~parfreeDistribution
    error('Distribution model (4th input) must either be a function handle, ''P'', or ''none''.')
end

% Get information about background parameters
% ===========================================
isparamodel = cellfun(@(model) isa(model,'function_handle'),bg_model);
includeBackground = ~cellfun(@(model) strcmp(model,'none'),bg_model);

[par0_bg,lower_bg,upper_bg,N_bg] = getmodelparams(bg_model);

if any(~isparamodel & includeBackground)
    error('Background model (5th input) must either be a function handle, or ''none''.')
end

    
% Get information about experiment parameters
% ===========================================
isparamodel = cellfun(@(model) isa(model,'function_handle'),ex_model);
includeExperiment = ~cellfun(@(model) strcmp(model,'none'),ex_model);

[par0_ex,lower_ex,upper_ex,N_ex] = getmodelparams(ex_model);

if  any(~isparamodel & includeExperiment)
    error('Experiment models must either be a function handle, or ''none''.')
end

% Catch nonsensical situation
if any(~includeForeground & ~includeExperiment)
    error('Cannot fit anything without distribution model and without background model.')
end

% Prepare full-parameter set
% ===========================
% Combine all initial, lower and upper bounds of parameters into vectors
par0 = parcombine(par0,{par0_dd,par0_bg,par0_ex});
lb = parcombine(lb,{lower_dd,lower_bg,lower_ex});
ub = parcombine(ub,{upper_dd,upper_bg,upper_ex});
checkbounds(lb,ub,par0);
Nparam = numel(par0);

% Build index vectors for accessing parameter subsets
ddidx = 1:N_dd;
[bgidx,exidx] = deal(cell(nSignals,1));
for ii = 1:nSignals
    bgidx{ii} = N_dd + sum(N_bg(1:ii-1)) + (1:N_bg(ii));
    exidx{ii} = N_dd + sum(N_bg) + sum(N_ex(1:ii-1)) + (1:N_ex(ii));
end


% Fit the dipolar multi-pathway model
% ====================================
OnlyRegularization = ~includeExperiment & ~includeBackground;
OnlyParametric = ~OnlyRegularization & (parametricDistribution || ~includeForeground);


if OnlyRegularization
    
    % Use basic dipolar kernel
    Ks = cellfun(@(t) dipolarkernel(t,r),t,'UniformOutput',false);
    
    % Linear regularization fit
    [Pfit,Pfit_uq] = fitregmodel(Vexp,Ks,r,regtype,regparam,'Verbose',verbose);
    
    % Get fitted models
    Vfit = cellfun(@(K) K*Pfit,Ks,'UniformOutput',false);
    Bfit = cellfun(@(V) ones(size(V)),Vexp,'UniformOutput',false);
    
    % No parameters
    parfit_ = [];
    [Vfit_uq,Bfit_uq,paruq_bg,paruq_dd,paruq_ex] = deal([]);
    
elseif OnlyParametric
    
    % Prepare the full-parametric model
    if includeForeground
        Pfcn = @(par) dd_model(r,par(ddidx));
    else
        Pfcn = @(~) ones(numel(r),1);
    end
    Vmodel = @(~,par) cellfun(@(K) K*Pfcn(par),multiPathwayModel(par),'UniformOutput',false);
    
    % Non-linear parametric fit
    [parfit_,Vfit,param_uq] = fitparamodel(Vexp,Vmodel,t,par0,lb,ub);

    % Get fitted models
    [~,Bfit] = multiPathwayModel(parfit_);
    if includeForeground
        Pfit = Pfcn(parfit_);
    else
        Pfit = [];
    end
    if ~iscell(Vfit)
        Vfit = {Vfit};
    end
    if calculateCI
    [Vfit_uq,Pfit_uq,Bfit_uq,paruq_bg,paruq_ex,paruq_dd] = splituq(param_uq);
    end
    
else
    
    % Non-negativity constraint on distributions
    lbl = zeros(numel(r),1);
    
    % Separable non-linear least squares (SNNLS) 
    [parfit_,Pfit,snlls_uq] = snlls(Vexp,@multiPathwayModel,par0,lb,ub,lbl);

    % Get the fitted models
    [Kfit,Bfit] = multiPathwayModel(parfit_);
    Pfit = Pfit(:);
    Vfit = cellfun(@(K) K*Pfit,Kfit,'UniformOutput',false);
    
    if calculateCI
    [Vfit_uq,Pfit_uq,Bfit_uq,paruq_bg,paruq_ex,paruq_dd] = splituq(snlls_uq);
    end
end

% Normalize distribution
% =======================
if normP && includeForeground
    Pnorm = trapz(r,Pfit);
    Pfit = Pfit/Pnorm;
    if calculateCI
        % scale CIs accordingly
        if iscell(Pfit_uq)
            for j = 1:numel(Pfit_uq)
                Pfit_uq{j} = Pfit_uq{j}/Pnorm;
            end
        else
            Pfit_uq.ci = @(p) Pfit_uq.ci(p)/Pnorm;
        end
    end
end

% Calculate goodness of fit
% =========================
if computeStats
    stats = cell(nSignals,1);
    for j = 1:nSignals
        Ndof = numel(Vexp{j}) - numel(parfit_);
        stats{j} = gof(Vexp{j},Vfit{j},Ndof);
    end
else
    stats = [];
end

% Return fitted parameters and confidence intervals in structures
% ================================================================
parfit_ = parfit_(:);
parfit.dd = parfit_(ddidx);
for j = 1:nSignals
    parfit.bg{j} = parfit_(bgidx{j});
    parfit.ex{j} = parfit_(exidx{j});
end
if calculateCI
    paruq.dd = paruq_dd;
    modfituq.Pfit = Pfit_uq;
    for j = 1:nSignals
        paruq.bg{j} = paruq_bg{j};
        paruq.ex{j} = paruq_ex{j};
        modfituq.Vfit{j} = Vfit_uq;
        modfituq.Bfit{j} = Bfit_uq;
    end
end

% Plotting
% =========
if DisplayResults
    display()
end

% Return numeric and not cell arrays if there is only one signal
if nSignals==1
    if iscell(Vfit)
    Vfit = Vfit{1};
    end
    if iscell(Bfit)
    Bfit = Bfit{1};
    end
    if iscell(parfit.dd)
        parfit.dd = parfit.dd{1};
    end
    if iscell(parfit.bg)
        parfit.bg = parfit.bg{1};
    end
    if iscell(parfit.ex)
        parfit.ex = parfit.ex{1};
    end
    if calculateCI
        paruq.bg = paruq.bg{1};
        paruq.ex = paruq.ex{1};
    end
    if ~isempty(stats)
        stats = stats{1};
    end
end

% =========================================================================
% Subfunctions
% =========================================================================

% Multi-pathway dipolar model
% ------------------------------------------------------------------
% This function represent the core model of fitsignal, it takes the
% parameters and computes the dipolar kernel and background according to
% the dipolar multi-pathway theory. 
% This is the non-linear part of the SNLLS problem when fitting a
% parameter-free distribution.
    function [Ks,Bs] = multiPathwayModel(par)
        
        [Bs,Ks] = deal(cell(nSignals,1));
        for iSignal=1:nSignals            
            % Get parameter subsets for this signal
            ex_par = par(exidx{iSignal});
            bg_par = par(bgidx{iSignal});
            
            % Prepared background basis function
            if includeBackground(iSignal)
                Bfcn = @(t,lam) bg_model{iSignal}(t,bg_par,lam);
            else
                Bfcn = @(~,~) ones(size(Vexp{iSignal}));
            end
            
            % Get pathway information
            if includeExperiment(iSignal)
                pathinfo = ex_model{iSignal}(ex_par);
            else
                pathinfo = 1;
            end
            
            % Compute the multipathway-background
            Bs{iSignal} = dipolarbackground(t{iSignal},pathinfo,Bfcn);
            % Compute the multipathway-kernel
            Ks{iSignal} = dipolarkernel(t{iSignal},r,pathinfo);
            Ks{iSignal} = Ks{iSignal}.*Bs{iSignal};
            
        end
        
    end
    function [Ks] = multiPathwayKernel(par,idx)
        Ks = multiPathwayModel(par);
        if nargin>1
            Ks = Ks{idx};
        end
    end


    function [Bs] = multiPathwayBackground(par,idx)
        [~,Bs] = multiPathwayModel(par);
        if nargin>1
            Bs = Bs{idx};
        end
    end

% Uncertainty quantification
% ------------------------------------------------------------------
% This function takes the full combined covariance matrix of the
% linear+nonlinear parameter sets and splits it into the different
% components of the model. 
    function [Vfit_uq,Pfit_uq,Bfit_uq,paruq_bg,paruq_ex,paruq_dd] = splituq(full_uq)
        
        % Pre-allocation
        [paruq_bg,paruq_ex,Bfit_uq,Vfit_uq] = deal(cell(nSignals,1));
        
        % Retrieve full covariance matrix
        covmat = full_uq.covmat;
        
        paramidx = 1:Nparam;
        Pfreeidx = Nparam+1:length(covmat);
        
        % Full parameter set uncertainty
        % ==================================
        subcovmat = covmat(paramidx,paramidx);
        paruq = uqst('covariance',parfit_,subcovmat,lb,ub);
        
        % Background parameters uncertainty
        % ==================================
        for jj=1:nSignals
            if includeBackground(jj)
                bgsubcovmat  = paruq.covmat(bgidx{jj},bgidx{jj});
                paruq_bg{jj} = uqst('covariance',parfit_(bgidx{jj}),bgsubcovmat,lb(bgidx{jj}),ub(bgidx{jj}));
            else
                paruq_bg{jj} = [];
            end
        end
        
        % Experiment parameters uncertainty
        % ==================================
        for jj=1:nSignals
            if includeExperiment(jj)
                exsubcovmat  = paruq.covmat(exidx{jj},exidx{jj});
                paruq_ex{jj} = uqst('covariance',parfit_(exidx{jj}),exsubcovmat,lb(exidx{jj}),ub(exidx{jj}));
            else
                paruq_ex{jj}  = [];
            end
        end
        
        % Distribution parameters uncertainty
        % ====================================
        if parametricDistribution
            ddsubcovmat  = paruq.covmat(ddidx,ddidx);
            paruq_dd = uqst('covariance',parfit_(ddidx),ddsubcovmat,lb(ddidx),ub(ddidx));
        else
            paruq_dd = [];
        end
        
        % Distance distribution uncertainty
        % ====================================
        nonneg = zeros(numel(r),1);
        if parametricDistribution
            Pfit_uq = paruq.propagate(Pfcn,nonneg,[]);
        else
            subcovmat = covmat(Pfreeidx,Pfreeidx);
            Pfit_uq = uqst('covariance',Pfit,subcovmat,nonneg,[]);
        end
        
        % Background uncertainty
        % ====================================
        for jj=1:nSignals
            if includeExperiment(jj)
                Bfit_uq{ii} = paruq.propagate(@(par)multiPathwayBackground(par,jj));
            else
                Bfit_uq{ii} = [];
            end
        end
        
        % Dipolar signal uncertainty
        % ====================================
        for jj=1:nSignals
            if parametricDistribution
                % Simple parametric model error propagation
                Vmodel = @(par)multiPathwayKernel(par,jj)*Pfcn(par(ddidx));
                Vfit_uq{ii} = paruq.propagate(Vmodel);
            else
                % Use special structure to speed up propagation for
                % parameter-free case instead of .propagate()
                J = [jacobianest(@(par)multiPathwayKernel(par(paramidx),jj)*Pfit,parfit_) Kfit{jj}];
                Vcovmat = J*covmat*J.';
                Vfit_uq{ii} = uqst('covariance',Vfit{jj},Vcovmat,[],[]);
            end
        end
        
    end

% Parsing and validation of options
% ------------------------------------------------------------------
% This function parses the name-value pairs or structures with options
% passes to the main function in the varargin. If the options pass, the
% default values are overwritten by the user-specified values.
    function parsevalidate(varargin)
        
        validoptions = {'regparam','regtype','normP','Display'};
        
        % Parse options
        [regparam_,regtype_,normP_,DisplayResults_] = parseoptions(validoptions,varargin);

        if ~isempty(DisplayResults_)
            validateattributes(DisplayResults_,{'logical'},{'scalar','nonempty'})
            DisplayResults = DisplayResults_;
        end
        if ~isempty(normP_)
            validateattributes(normP_,{'logical'},{'scalar','nonempty'})
            normP = normP_;
        end
        if ~isempty(regparam_)
            if ischar(regparam)
                allowedMethodInputs = {'lr','lc','cv','gcv','rgcv','srgcv','aic','bic','aicc','rm','ee','ncp','gml','mcl'};
                validateattributes(regparam,{'char'},{'nonempty'},mfilename,'RegParam option');
                regparam = validatestring(regparam,allowedMethodInputs);
            else
                validateattributes(regparam,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'RegParam')
            end
        end
        if ~isempty(regtype_)
            validateattributes(regtype_,{'char'},{'scalar','nonempty'})
            regtype = validatestring(regtype_,{'tikh','tv','huber'});
        end
    end

% Validation of mandatory inputs
% ------------------------------------------------------------------
% This function validates the inputs signals, time-axes and distance axis
    function nSignals = validateVtr()
        if ~iscell(Vexp)
            Vexp = {Vexp};
        end
        nSignals = numel(Vexp);
        r = r(:).';
        if ~iscell(t)
            t = {t};
        end
        if numel(t)~=nSignals
            error('The same number of signals V and time axes t must be provided.')
        end
        for i = 1:nSignals
            Vexp{i} = Vexp{i}(:);
            t{i} = t{i}(:);
            if numel(Vexp{i})~=numel(t{i})
                error('V{%i} and t{%i} must have the same number of elements.',i,i)
            end
            validateattributes(Vexp{i},{'numeric'},{'vector'},mfilename,'V (1st input)');
            validateattributes(t{i},{'numeric'},{'vector'},mfilename,'t (2nd input)');
            if ~isreal(Vexp{i})
                error('Input signals cannot be complex-valued')
            end
        end
        validateattributes(r,{'numeric'},{'vector'},mfilename,'r (3rd input)');
        
        if ~iscell(bg_model)
            bg_model = {bg_model};
        end
        if numel(bg_model)~=nSignals
            bg_model = repmat(bg_model,nSignals,1);
        end
        
        
        if ~iscell(ex_model)
            ex_model = {ex_model};
        end
        if numel(ex_model)~=nSignals
            ex_model = repmat(ex_model,nSignals,1);
        end
        
    end

    function display()
        for i = 1:nSignals
            subplot(2,nSignals,i);
            plot(t{i},Vexp{i},'k.',t{i},Vfit{i},'r','LineWidth',1.5)
            hold on
            Vci95 = Vfit_uq{i}.ci(95);
            Vci50 = Vfit_uq{i}.ci(50);
            fill([t{i}; flipud(t{i})],[Vci95(:,1); flipud(Vci95(:,2))],'r','LineStyle','none','FaceAlpha',0.2)
            fill([t{i}; flipud(t{i})],[Vci50(:,1); flipud(Vci50(:,2))],'r','LineStyle','none','FaceAlpha',0.5)
            hold off
            axis tight
            grid on
            xlabel('Time [\mus]');
            ylabel(sprintf('V\\{%d\\}',i));
            legend('Data','Fit','95%-CI','50%-CI');
        end
        subplot(2,1,2);
        plot(r,Pfit,'k','LineWidth',1.5);
        hold on
        Pci95 = Pfit_uq.ci(95);
        Pci50 = Pfit_uq.ci(50);
        fill([r fliplr(r)],[Pci95(:,1); flipud(Pci95(:,2))],'r','LineStyle','none','FaceAlpha',0.2)
        fill([r fliplr(r)],[Pci50(:,1); flipud(Pci50(:,2))],'r','LineStyle','none','FaceAlpha',0.5)
        hold off
        xlabel('Distance [nm]');
        ylabel('P [nm^{-1}]');
        legend('Fit','95%-CI','50%-CI')
        axis tight,grid on
        drawnow
        disp('Goodness of fit')
        for i=1:nSignals
            fprintf('  V{%i}: %s2 = %.4f  RMSD  = %.4f \n',i,char(hex2dec('3c7')),stats{i}.chi2red,stats{i}.RMSD)
        end
        disp('Fitted parameters and 95%-confidence intervals')
        str = '  %s{%d}(%d):   %.7f  (%.7f, %.7f)  %s (%s)\n';
        if numel(parfit.dd)>0
            info = dd_model();
            ci = paruq.dd.ci(95);
            for p = 1:numel(parfit.dd)
                c = parfit.dd(p);
                fprintf(str,'dd',1,p,c,...
                    ci(p,1),ci(p,2),info(p).Parameter,info(p).Units)
            end
        end
        if numel(parfit.bg)>0
            for i = 1:nSignals
                info = bg_model{i}();
                ci = paruq.bg{i}.ci(95);
                for p = 1:numel(parfit.bg{i})
                    c = parfit.bg{i}(p);
                    fprintf(str,'bg',i,p,c,...
                        ci(p,1),ci(p,2),info(p).Parameter,info(p).Units)
                end
            end
        end
        if numel(parfit.ex)>0
            for i = 1:nSignals
                info = ex_model{i}();
                ci = paruq.ex{i}.ci(95);
                for p = 1:numel(parfit.ex{i})
                    c = parfit.ex{i}(p);
                    fprintf(str,'ex',i,p,c,...
                        ci(p,1),ci(p,2),info(p).Parameter,info(p).Units)
                end
            end
        end
    end

end

%===============================================================================
% Local functions
%===============================================================================

% Extract parameter info from parametric models
% ------------------------------------------------------------------
function [par0,lo,up,N] = getmodelparams(models)

if ~iscell(models)
    models = {models};
end

[par0,lo,up] = deal(cell(numel(models),1));
N = zeros(numel(models),1);
for i=1:numel(models)
    if isa(models{i},'function_handle')
        info = models{i}();
        par0{i} = [info.Start];
        lo{i} = [info.Lower];
        up{i} = [info.Upper];
    else
        [par0{i},lo{i},up{i}] = deal([]);
    end
    N(i) = numel(par0{i});
end

if numel(models)==1
    par0 = par0{1};
    lo = lo{1};
    up = up{1};
end

end

% Combine parameter subsets
% ------------------------------------------------------------------
function pvec = parcombine(p,def)
for k = 1:3
    if isempty(p{k}), p{k} = def{k}; end
    if iscell(p{k}), p{k} = cell2mat(p{k}); end
    p{k} = p{k}(:).';
end
pvec = cell2mat(p);
end


% Boundary checking
% ------------------------------------------------------------------
function checkbounds(lb,ub,par0)

nParams = numel(par0);
if numel(lb)~=nParams
    error('Lower bounds and initial values of parameters must have the same number of elements.');
end
if numel(ub)~=nParams
    error('Lower bounds and initial values of parameters must have the same number of elements.');
end
if any(lb>ub)
    error('Lower bounds cannot be larger than upper bounds.');
end
if any(lb>par0 | par0>ub)
    error('Inital values for parameters must lie between lower and upper bounds.');
end

end

% Input parsing
% ------------------------------------------------------------------
% This function parses the different input schemes
function [dd_model,bg_model,ex_model,par0,lb,ub,options] = parseinputs(inputs)

% Prepare empty containers
[dd_model,bg_model,ex_model,par0,lb,ub] = deal([]);
% Parse the ipnuts cell array
optionstart = numel(inputs);
for i=1:numel(inputs)
    if i<4 && ischar(inputs{i}) && ~any(strcmp(inputs{i},{'P','none'})) || ...
            i>3 && ischar(inputs{i})
        optionstart = i-1;
        break
    end
    switch i
        case 1
            dd_model = inputs{1};
        case 2
            bg_model = inputs{2};
        case 3
            ex_model = inputs{3};
        case 4
            par0 = inputs{4};
        case 5
            lb = inputs{5};
        case 6
            ub = inputs{6};
    end
end
inputs(1:optionstart) = [];
options = inputs;

% Set defaults
if isempty(dd_model), dd_model = 'P'; end
if isempty(bg_model), bg_model = @bg_hom3d; end
if isempty(ex_model), ex_model = @ex_4pdeer; end

%Validate inputs
if isempty(lb)
    lb = {[],[],[]};
end
validateattributes(lb,{'cell'},{'row'},mfilename,'Lower option');
if isempty(ub)
    ub = {[],[],[]};
end
validateattributes(ub,{'cell'},{'row'},mfilename,'Upper option');

if isempty(par0)
    par0 = {[],[],[]};
else
    if ~iscell(par0) || numel(par0)~=3
        error('Initial parameters (7th input) must be a 3-element cell array.')
    end
end

end




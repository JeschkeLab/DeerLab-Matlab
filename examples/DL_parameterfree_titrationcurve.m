%============================================================================
% DeerLab Example:
% Analyzing pseudo-titration (dose-respononse) curves with 
% parameter-free distributions 
%============================================================================

% This example shows how to use separable non-linear least squares (SNLLS)
% to fit a pseudo-titration curve to multiple DEER datsets, using
% parameter-free distance distributions.

clc,clear,clf

% Generating multiple datasets
%-----------------------------------------------------------------------------

% First, let's start by simulating some data. In this example we will
% simulate a protein system in their states A (natural) and B (changed upon addition
% of a ligand L) given by the chemical equilibrium  A + L <-> B.

% Time axes
ts{1} = linspace(-0.2,3,100);
ts{2} = linspace(-0.1,5,300);
ts{3} = linspace(-0.5,2,200);
ts{4} = linspace(-0.1,1,100);
ts{5} = linspace(-0.2,6,300);
ts{6} = linspace(-0.2,3,300);
ts{7} = linspace(-0.1,4,100);
Nsignals = numel(ts);

% Distance axes for states A and B
rA = linspace(1,8,100);
rB = linspace(1,8,100);

% Distributions for states A and B
PA = dd_gauss(rA,[5.5 0.4]);
PB = dd_gauss2(rB,[4.5 0.7 0.4 3.5 0.6 0.6]);

L = [0.3 1 3 10 30 100 300]; % total ligand concentration, uM
KD = 5.65;  % dissociation constant, uM

% Populations of states A and B
[xA,xB] = chemicalequilibrium(KD,L);

% Global kernel model
Ks = Kmodel([0.25 0.1 KD],ts,rA,rB,L);

% Simulate dipolar signals
for i=1:Nsignals
    Vs{i} = Ks{i}*[PA; PB] + whitegaussnoise(ts{i},0.005);
    subplot(221)
    hold on,plot(ts{i},Vs{i}+i/9,'k.'),hold off
end

% Analysis
%-----------------------------------------------------------------------------

% For simplification, we will assume that all DEER traces have the same
% background function and modulation depth. Thus, we will fit the
% modulations depth (lam) and background decay constant (k) globally along
% the dissociation constant (KD).

% Non-linear parameters:
%       lam  k   KD
par0 = [0.5 0.5  5];  % start values 
lb   = [ 0   0   1];  % lower bounds
ub   = [ 1   1   10]; % upper bounds

% Linear parameters:
%     |-------PA--------||--------PB--------|
lbl = [zeros(1,numel(rA)) zeros(1,numel(rB))]; % Non-negativity constraint
ubl = []; % Unconstrained

% Run SNLLS optimization
[parfit,Pfit,puq] = snlls(Vs,@(p)Kmodel(p,ts,rA,rB,L),par0,lb,ub,lbl,ubl,'GlobalWeights',ones(Nsignals,1));

% Extract the fitted disociation constant value and its 95%-confidence interval
KD = parfit(3);
parci = puq.ci(95,'nonlin');
KDci = parci(3,:);

% Print result
fprintf('KD = %.2f(%.2f-%.2f)uM \n',KD,KDci(1),KDci(2))
% Plot results
plotresults(parfit,Pfit,ts,rA,rB,L,xA,xB)


% Dipolar kernel model
%-----------------------------------------------------------------------------
% Here we define the dipolar kernel model as the non-linear function of the
% SNLLS problem. This function needs to take the parameters and return a
% cell-array of kernels, each one for the corresponding datasets that we
% have. 
% Since we have a total distribution of the form 
%     P = xA*PA + xB*PB
% we can define an augmented kernel as
%     K = [xA*KA xB*KB]
% such that 
%     K*[PA PB] = V
% and the vector [PA PB] constitutes the linear part fitted by SNLLS.
function K = Kmodel(p,ts,rA,rB,L)

Nsignals = numel(ts);

% Unpack parameters
lam = p(1); 
k   = p(2);
KD = p(3);

% Get fractions for given KD
[xA,xB] = chemicalequilibrium(KD,L);

% General the dipolar kernels
for i=1:Nsignals
    B = bg_exp(ts{i},k,lam);
    % Kernel for fraction A
    KA = dipolarkernel(ts{i},rA,lam,B);
    % Kernel for fraction B
    KB = dipolarkernel(ts{i},rB,lam,B);
    K{i} = [xA(i)*KA xB(i)*KB];
end

end


%Prepare equilibrium of type: A + L <-> B
%----------------------------------------------------------------------------

function [xA,xB] = chemicalequilibrium(KD,L)

Ctot = 1; % total protein concentration, uM

% % Get fraction of state B
Kb = 1/KD;
for q = 1:numel(L)
    xB_ = roots([Kb*Ctot -(Kb*L(q) + Kb*Ctot + 1) Kb*L(q)]);
    try
    xB(q) = xB_(xB_<=1 & xB_>=0);
    catch
    xB(q) = min(1,max(0,xB_(1)));    
    end
end
% Get fraction of state A
xA = 1 - xB;
end


% Plots
%-----------------------------------------------------------------------------

function plotresults(parfit,Pfit,ts,rA,rB,L,xA,xB)

Nsignals = numel(ts);

% Simulate fits
Ksfit = Kmodel(parfit,ts,rA,rB,L);
for i=1:Nsignals
    Vsfit{i} = Ksfit{i}*Pfit(:);
    subplot(221)
    hold on,plot(ts{i},Vsfit{i}+i/9,'r'),hold off
end
axis tight, grid on, box on
xlabel('t [\mus]')
ylabel('V(t)')
legend('data','fit')
KDfit = parfit(3);
[xAfit,xBfit] = chemicalequilibrium(KDfit,L);


for i=1:Nsignals
    PAfit = xAfit(i)*Pfit(1:numel(rA));
    PBfit = xBfit(i)*Pfit(numel(rA)+1:end);
    subplot(2,2,[2,4])
    hold on
    plot(rA,PAfit+2*i,'r',rB,PBfit+2*i,'b')
end
axis tight, grid on, box on
xlabel('r [nm]')
ylabel('P(r)')
legend('state A','state B')

subplot(223)
cla,hold on
plot(log10(L),xA,'r-',log10(L),xB,'b-')
plot(log10(L),xAfit,'ro',log10(L),xBfit,'bo')
axis tight, grid on, box on
xlabel('log_{10}([L])')
ylabel('Fractions')
legend('state A','state B')
end



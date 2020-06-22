% =================================================================================
%  MINQ5
% =================================================================================
%  W. Huyer and A. Neumaier
%  MINQ8 - General Definite and Bound Constrained Indefinite Quadratic Programming
%  Manuscript (2017). http://www.mat.univie.ac.at/~neum/software/minq/ 
% =================================================================================
%
% Downloaded from https://www.mat.univie.ac.at/~neum/software/minq/ on 06/22/2020
%

% CONTENT OF MINQ'S LICENSE
% -------------------------------------------------------------------------------
% THIS TRANSLATION IS FOR YOUR INFORMATION ONLY!
% THE LEGALLY BINDING DOCUMENT IS THE ORIGINAL LICENSE
% AGREEMENT IN GERMAN,
% Please find the German version of the license agreement at the end of this
% document.
% License Agreement on the MINQ Software provided by Prof. Arnold Neumair.
% Version 1.0 – December 2016
% 1. Objective and scope of the License Agreement
% 2. Usage by end user
% 3. Distribution and modification
% 4. Revocation of rights of use
% 5. Liability of the right holder
% 6. Final provisions
% 1. Objective and scope of the License Agreement
% 1.1
% The MINQ Software Package is made available to you free of charge by Prof.
% Arnold Neumaier, and may be used for private as well as for academic purposes.
% 1.2 The user has no right to claim the right to use of the MINQ Software Package
% or to its provision. The software and - as the case may be – any accompanying
% documentation are provided “as is”. Prof. Arnold Neumaier does not guarantee for
% the proper function of the Software or its fitness for a particular purpose and does
% not provide support. In this respect only the resources offered at
% http://www.mat.univie.ac.at/~neum/software/minq/ are available.
% 1.3 Also in the future the existing functions of the MINQ shall remain free of
% charge. However, the owner reserves the right to integrate additional components
% to the program which might be subject to fees.
% 2. Usage by end user
% 2.1 The MINQ can be used by the end user (private or academic) for an
% unspecified time.
% 2.2 For commercial use of the MINQ Software please contact Prof. Arnold
% Neumaier under (Arnold.Neumaier@univie.ac.at).
% 3. Distribution and modification
% 3.1 The right to claim damages is reserved.
% 3.2 The program may only be distributed free of charge and only together with all
% accompanying files from the site
% http://www.mat.univie.ac.at/~neum/software/minq/
% 3.3 Any use of the MINQ in a way or in a context which might cause damage to
% the reputation of the project MINQ, the University of Vienna, or to Prof. Arnold
% Neumaier is prohibited.
% 4. Revocation of rights of use
% The right holder has the right to revoke all rights of use of the user in case the
% user is in breach of this License Agreement. The user admits that all used product
% names and registered brands / trade mark are the property of their respective
% owners, no matter whether they are marked as those or not.
% 5. Liability of the right holder
% 5.1 Prof. Arnold Neumaier is not liable for damages which were caused by the
% software, irrespective of the legal basis – neither of contractual nor of extracontractual
% nature. This restriction of liability also applies when services are
% offered for remuneration.
% 6. Final provisions
% 6.1 This License Agreement remains in force, until it is replaced by a newer
% version. Prof. Arnold Neumaier has the right to make the installation of a new
% program version depend on the users consent with a renewed License
% Agreement.
% 6.2 The legally binding version of this License Agreement is solely the German
% version. Only this German version is decisive for the content of this License
% Agreement and the rights and duties arising from it. Versions in other languages
% are non-binding translations which are merely for information purposes.
% 6.3 This License Agreement is governed by the laws of the Republic of Austria to
% the exclusion of the Convention on the International Sale of Goods. The place of
% jurisdiction is Vienna. In any case, the parties can be sued in their respective
% place of general jurisdiction.
%
% -------------------------------------------------------------------------------
%
% Lizenzbedingungen für MINQ von Prof. Arnold Neumaier
% Version 1.0 – Dezember 2016
% 1. Gegenstand der Lizenzvereinbarung
% 2. Nutzung durch Endanwender
% 3. Weitergabe und Veränderung
% 4. Entzug der Nutzungsrechte
% 5. Haftung des Rechteinhabers
% 6. Schlussbestimmungen
% 1. Gegenstand der Lizenzvereinbarung
% 1.1 MINQ wird ihnen kostenlos von Prof. Arnold Neumaier, zur Verfügung gestellt
% und darf sowohl für private als auch für universitäre Zwecke genutzt werden.
% 1.2 Der Nutzer hat keinen Anspruch auf die Nutzung des MINQ Pakets oder
% dessen Zurverfügungstellung. Die Software und eine ggf. zugehörige
% Dokumentation werden so wie sie sind zur Verfügung gestellt. Prof. Arnold
% Neumaier übernimmt keine Gewähr für die ordnungsgemäße Funktion des MINQ
% Paketes oder die Erreichung eines bestimmten Zwecks und leistet keinen
% Support. Insoweit stehen lediglich die unter
% http://www.mat.univie.ac.at/~neum/software/minq/ angebotenen Ressourcen zur
% Verfügung.
% 1.3 Die bestehenden Funktionen von MINQ sollen auch in Zukunft kostenfrei
% bleiben. Prof. Arnold Neumaier behält sich jedoch die Integration zusätzlicher
% kostenpflichtiger Programmteile vor.
% 2. Nutzung durch Endanwender
% 2.1 MINQ steht dem privaten und universitären Anwender ab Installation für unbestimmte
% Zeit, bis auf Wiederruf, zur Verfügung.
% 2.2 Für gewerbliche Anwendungen der MINQ Software kontaktieren Sie bitte
% Prof. Arnold Neumaier unter (Arnold.Neumaier@univie.ac.at).
% 3. Veränderung und Weitergabe
% 3.1 Die Geltendmachung von Schadensersatz wird ausdrücklich vorbehalten.
% 3.2 Das Programm darf nur kostenlos und nur zusammen mit allen zugehörigen
% Dateien und in unverändertem Zustand auf
% http://www.mat.univie.ac.at/~neum/software/minq/ verbreitet werden.
% 3.3 Eine Weiterverbreitung von MINQ in einer Art und Weise oder in einem Zusammenhang,
% die dem Ansehen des Projekts MINQ, oder der Universität Wien,
% deren Angestellten, Vertretern, Erfüllungsgehilfen oder Partnern und Prof. Arnold
% Neumaier schaden kann, ist unzulässig.
% 4. Entzug der Nutzungsrechte
% Bei Verstoß gegen diese Lizenzbedingungen ist der Rechteinhaber berechtigt,
% dem Nutzer alle Nutzungsrechte zu entziehen. Alle verwendeten Produktnamen
% und eingetragenen Marken/Warenzeichen werden hiermit als Eigentum ihrer
% Inhaber anerkannt, unabhängig davon, ob sie als solche gekennzeichnet sind
% oder nicht.
% 5. Haftung des Rechteinhabers
% 5.1 Prof. Arnold Neumaier haftet nicht für Schäden, die durch die angebotene
% Software verursacht wurde, gleich aus welchem Rechtsgrund - sowohl
% vertraglicher als auch außervertraglicher. Diese Haftungsbeschränkung gilt auch
% bei gegen Entgelt erbrachten Diensten.
% 6. Schlussbestimmungen
% 6.1 Diese Lizenzbedingungen bleiben gültig, bis sie durch eine neuere Version
% ersetzt werden. Prof. Arnold Neumaier ist berechtigt, die Installation einer neuen
% Programmversion von der Zustimmung des Anwenders zu geänderten
% Lizenzbedingungen abhängig zu machen.
% 6.2 Verbindliche Fassungen dieser Lizenzbedingungen sind ausschließlich die
% Fassung in deutscher Sprache. Nur diese Fassung ist für den Inhalt dieser
% Lizenzbedingungen und die Rechte und Pflichten aus ihnen maßgeblich.
% Fassungen in anderen Sprachen sind unverbindliche Übersetzungen, die lediglich
% Informationszwecken dienen.
% 6.3 Es gilt das Recht der Republik Österreich unter Ausschluss des UN-Kaufrechts.
% Gerichtsstand ist Wien. Die Parteien können jeweils auch an ihrem allgemeinen
% Gerichtsstand verklagt werden.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% minq.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,fct,ier,nsub]=minq(gam,c,G,xu,xo,prt,xx);
% minimizes an affine quadratic form subject to simple bounds,
% using coordinate searches and reduced subspace minimizations
% using LDL^T factorization updates
%    min    fct = gam + c^T x + 0.5 x^T G x
%    s.t.   x in [xu,xo]    % xu<=xo is assumed
% where G is a symmetric n x n matrix, not necessarily definite
% (if G is indefinite, only a local minimum is found)
%
% if G is sparse, it is assumed that the ordering is such that
% a sparse modified Cholesky factorization is feasible
%
% prt	printlevel
% xx	initial guess (optional)
%
% x	minimizer (but unbounded direction if ier=1)
% fct	optimal function value
% ier	0  (local minimizer found)
% 	1  (unbounded below)
% 	99 (maxit exceeded)
% 	-1 (input error)
%
% calls getalp.m, ldl*.m, minqsub.m, pr01.m
%
function [x,fct,ier,nsub]=minq(gam,c,G,xu,xo,prt,xx)




%%% initialization %%%
%disp('####################################################');
%disp('######################  minq  ######################');
%disp('####################################################');
convex=0;
n=size(G,1);

if prt>0,
    % show printlevel
    printlevel=prt
    
    % check input data for consistency
    ier=0;
    if size(G,2)~=n,
        ier=-1;disp('minq: Hessian has wrong dimension');
        x=NaN+zeros(n,1);fct=NaN;nsub=-1;
        return;
    end;
    if size(c,1)~=n | size(c,2)~=1,
        ier=-1;disp('minq: linear term has wrong dimension');
    end;
    if size(xu,1)~=n | size(xu,2)~=1,
        ier=-1;disp('minq: lower bound has wrong dimension');
    end;
    if size(xo,1)~=n | size(xo,2)~=1,
        ier=-1;disp('minq: lower bound has wrong dimension');
    end;
    if exist('xx')==1,
        if size(xx,1)~=n | size(xx,2)~=1,
            ier=-1;disp('minq: lower bound has wrong dimension');
        end;
    end;
    if ier==-1,
        x=NaN+zeros(n,1);fct=NaN;nsub=-1;
        return;
    end;
end;

maxit=3*n;       	% maximal number of iterations
% this limits the work to about 1+4*maxit/n matrix multiplies
% usually at most 2*n iterations are needed for convergence
nitrefmax=3;		% maximal number of iterative refinement steps

% initialize trial point xx, function value fct and gradient g

if nargin<7,
    % cold start with absolutely smallest feasible point
    xx=zeros(n,1);
end;
% force starting point into the box
xx=max(xu,min(xx,xo));

% regularization for low rank problems
hpeps=100*eps;	% perturbation in last two digits
G=G+spdiags(hpeps*diag(G),0,n,n);

% initialize LDL^T factorization of G_KK
K=logical(zeros(n,1));	% initially no rows in factorization
if issparse(G), L=speye(n); else L=eye(n); end;
dd=ones(n,1);

% dummy initialization of indicator of free variables
% will become correct after first coordinate search
free=logical(zeros(n,1));
nfree=0;
nfree_old=-1;

fct=inf; 		% best function value
nsub=0;			% number of subspace steps
unfix=1;		% allow variables to be freed in csearch?
nitref=0;		% no iterative refinement steps so far
improvement=1;		% improvement expected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop: alternating coordinate and subspace searches
while 1,
    if prt>1, disp('enter main loop'); end;
    if norm(xx,inf)==inf, error('infinite xx in minq.m'); end;
    g=G*xx+c;
    fctnew=gam+0.5*xx'*(c+g);
    if ~improvement,
        % good termination
        if prt,
            disp('terminate: no improvement in coordinate search');
        end;
        ier=0; break;
    elseif nitref>nitrefmax,
        % good termination
        if prt, disp('terminate: nitref>nitrefmax'); end;
        ier=0; break;
    elseif nitref>0 & nfree_old==nfree & fctnew >= fct,
        % good termination
        if prt,
            disp('terminate: nitref>0 & nfree_old==nfree & fctnew>=fct');
        end;
        ier=0; break;
    elseif nitref==0,
        x=xx;
        fct=min(fct,fctnew);
        if prt>1, fct, end;
        if prt>2, X=x', fct, end;
    else % more accurate g and hence f if nitref>0
        x=xx;
        fct=fctnew;
        if prt>1, fct, end;
        if prt>2, X=x', fct, end;
    end;
    if nitref==0 & nsub>=maxit,
        if prt,
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            disp('!!!!!           minq          !!!!!');
            disp('!!!!! incomplete minimization !!!!!');
            disp('!!!!!   too many iterations   !!!!!');
            disp('!!!!!     increase maxit      !!!!!');
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        else
            disp('iteration limit exceeded');
        end;
        ier=99;
        break;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % coordinate search
    count=0; 	% number of consecutive free steps
    k=0;     	% current coordinate searched
    while 1,
        while count<=n,
            % find next free index (or next index if unfix)
            count=count+1;
            if k==n, k=0; end;
            k=k+1;
            if free(k) | unfix, break; end;
        end;
        if count>n,
            % complete sweep performed without fixing a new active bound
            break;
        end;
        q=G(:,k);
        alpu=xu(k)-x(k); alpo=xo(k)-x(k); % bounds on step
        
        % find step size
        [alp,lba,uba,ier]=getalp(alpu,alpo,g(k),q(k));
        if ier,
            x=zeros(n,1);
            if lba, x(k)=-1; else x(k)=1; end;
            if prt,
                gTp=g(k),pTGp=q(k),quot=pTGp/norm(G(:),inf)
                disp('minq: function unbounded below in coordinate direction');
                disp('      unbounded direction returned');
                disp('      possibly caused by roundoff');
            end;
            if prt>1,
                disp('f(alp*x)=gam+gam1*alp+gam2*alp^2/2, where');
                gam1=c'*x
                gam2=x'*(G*x)
                ddd=diag(G);
                min_diag_G=min(ddd)
                max_diag_G=max(ddd)
            end;
            return;
        end;
        xnew=x(k)+alp;
        if prt & nitref>0,
            xnew,alp
        end;
        
        if lba | xnew<=xu(k),
            % lower bound active
            if prt>2, disp([num2str(k), ' at lower bound']); end;
            if alpu~=0,
                x(k)=xu(k);
                g=g+alpu*q;
                count=0;
            end;
            free(k)=0;
        elseif uba | xnew>=xo(k),
            % upper bound active
            if prt>2, disp([num2str(k), ' at upper bound']); end;
            if alpo~=0,
                x(k)=xo(k);
                g=g+alpo*q;
                count=0;
            end;
            free(k)=0;
        else
            % no bound active
            if prt>2, disp([num2str(k), ' free']); end;
            if alp~=0,
                if prt>1 & ~free(k),
                    unfixstep=[x(k),alp],
                end;
                x(k)=xnew;
                g=g+alp*q;
                free(k)=1;
            end;
        end;
        
    end;
    % end of coordinate search
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nfree=sum(free);
    if (unfix & nfree_old==nfree),
        % in exact arithmetic, we are already optimal
        % recompute gradient for iterative refinement
        g=G*x+c;
        nitref=nitref+1;
        if prt>0,
            disp('optimum found; iterative refinement tried');
        end;
    else
        nitref=0;
    end;
    nfree_old=nfree;
    gain_cs=fct-gam-0.5*x'*(c+g);
    improvement=(gain_cs>0 | ~unfix);
    
    if prt,
        % print (0,1) profile of free and return the number of nonzeros
        nfree=pr01('csrch ',free);
    end;
    if prt, gain_cs, end;
    if prt>2, X=x', end;
    
    % subspace search
    xx=x;
    if ~improvement | nitref>nitrefmax,
        % optimal point found - nothing done
    elseif nitref>nitrefmax,
        % enough refinement steps - nothing done
    elseif nfree==0,
        % no free variables - no subspace step taken
        if prt>0,
            disp('no free variables - no subspace step taken')
        end;
        unfix=1;
    else
        % take a subspace step
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% minqsub.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % patch for minq.m containing the subspace search
        %
        
        nsub=nsub+1;
        
        if prt>0,
            fct_cs=gam+0.5*x'*(c+g);
            format='*** nsub = %4.0f fct = %15.6e fct_cs = %15.6e';
            disp(sprintf(format,nsub,fct,fct_cs));
        end;
        
        % downdate factorization
        for j=find(free<K)', 	% list of newly active indices
            [L,dd]=ldldown(L,dd,j);
            K(j)=0;
            if prt>10,
                disp('downdate');
                fact_ind=find(K)'
            end;
        end;
        
        % update factorization or find indefinite search direchtion
        definite=1;
        for j=find(free>K)', 	% list of newly freed indices
            % later: speed up the following by passing K to ldlup.m!
            p=zeros(n,1);
            if n>1, p(K)=G(find(K),j); end;
            p(j)=G(j,j);
            [L,dd,p]=ldlup(L,dd,j,p);
            definite=isempty(p);
            if ~definite,
                if prt,disp('indefinite or illconditioned step');end;
                break;
            end;
            K(j)=1;
            if prt>10,
                disp('update');
                fact_ind=find(K)'
            end;
        end;
        
        if definite,
            % find reduced Newton direction
            p=zeros(n,1);
            p(K)=g(K);
            p=-L'\((L\p)./dd);
            if prt>10,
                disp('reduced Newton step');
                fact_ind=find(K)'
            end;
        end;
        
        if prt>2, input('continue with return>'); end;
        
        % set tiny entries to zero
        p=(x+p)-x;
        ind=find(p~=0);
        if isempty(ind);
            % zero direction
            if prt,disp('zero direction');end;
            unfix=1;
            return;
        end;
        
        % find range of step sizes
        pp=p(ind);
        oo=(xo(ind)-x(ind))./pp;
        uu=(xu(ind)-x(ind))./pp;
        alpu=max([oo(pp<0);uu(pp>0);-inf]);
        alpo=min([oo(pp>0);uu(pp<0);inf]);
        if alpo<=0 | alpu>=0,
            error('programming error: no alp');
        end;
        
        % find step size
        gTp=g'*p;
        agTp=abs(g)'*abs(p);
        if abs(gTp)<100*eps*agTp,
            % linear term consists of roundoff only
            gTp=0;
        end;
        pTGp=p'*(G*p);
        if convex, pTGp=max(0,pTGp); end;
        if ~definite & pTGp>0,
            if prt, disp(['tiny pTGp = ',num2str(pTGp),' set to zero']); end;
            pTGp=0;
        end;
        [alp,lba,uba,ier]=getalp(alpu,alpo,gTp,pTGp);
        if ier,
            x=zeros(n,1);
            if lba, x=-p; else x=p; end;
            if prt,
                qg=gTp/agTp
                qG=pTGp/(norm(p,1)^2*norm(G(:),inf))
                lam=eig(G);
                lam1=min(lam)/max(abs(lam))
                disp('minq: function unbounded below');
                disp('  unbounded subspace direction returned')
                disp('  possibly caused by roundoff');
                disp('  regularize G to avoid this!');
            end;
            if prt>1,
                disp('f(alp*x)=gam+gam1*alp+gam2*alp^2/2, where');
                gam1=c'*x
                rel1=gam1/(abs(c)'*abs(x))
                gam2=x'*(G*x)
                if convex, gam2=max(0,gam2); end;
                rel2=gam2/(abs(x)'*(abs(G)*abs(x)))
                ddd=diag(G);
                min_diag_G=min(ddd)
                max_diag_G=max(ddd)
            end;
            return;
        end;
        unfix=~(lba | uba);  % allow variables to be freed in csearch?
        
        % update of xx
        for k=1:length(ind),
            % avoid roundoff for active bounds
            ik=ind(k);
            if alp==uu(k),
                xx(ik)=xu(ik);
                free(ik)=0;
            elseif alp==oo(k),
                xx(ik)=xo(ik);
                free(ik)=0;
            else
                xx(ik)=xx(ik)+alp*p(ik);
            end;
            if abs(xx(ik))==inf,
                ik,alp,p(ik),
                error('infinite xx in minq.m');
            end;
            
        end;
        nfree=sum(free);
        subdone=1;
        
        
        if ier, return; end;
    end;
    
    if prt>0,
        % print (0,1) profile of free and return the number of nonzeros
        nfree=pr01('ssrch ',free);
        disp(' ');
        if unfix & sum(nfree)<n,
            disp('bounds may be freed in next csearch');
        end;
    end;
    
    
end;
% end of main loop
if prt>0,
    fct
    disp('################## end of minq ###################');
end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% getalp.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [alp,lba,uba,ier]=getalp(alpu,alpo,gTp,pTGp);
% get minimizer alp in [alpu,alpo] for a univariate quadratic
%	q(alp)=alp*gTp+0.5*alp^2*pTGp
% lba	lower bound active
% uba	upper bound active
%
% ier	 0 (finite minimizer)
%	 1 (unbounded minimum)
%
function [alp,lba,uba,ier]=getalp(alpu,alpo,gTp,pTGp);

lba=0;
uba=0;

% determine unboundedness
ier=0;
if alpu==-inf & ( pTGp<0 | (pTGp==0 & gTp>0) ),
    ier=1; lba=1;
end;
if alpo==inf & (pTGp<0 | (pTGp==0 & gTp<0) ),
    ier=1; uba=1;
end;
if ier, alp=NaN; return; end;

% determine activity
if pTGp==0 & gTp==0,
    alp=0;
elseif pTGp<=0,
    % concave case minimal at a bound
    if alpu==-inf,     lba=0;
    elseif alpo== inf, lba=1;
    else               lba = (2*gTp+(alpu+alpo)*pTGp>0);
    end;
    uba = ~lba;
else
    alp=-gTp/pTGp;          % unconstrained optimal step
    lba = (alp <= alpu);    % lower bound active
    uba = (alp >= alpo);    % upper bound active
end;

if lba, alp=alpu; end;
if uba, alp=alpo; end;


% print?
if abs(alp)==inf, gTp,pTGp,alpu,alpo,alp,lba,uba,ier, end;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ldlup.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [L,d,p]=ldlup(L,d,j,g)
% updates LDL^T factorization when a unit j-th row and column
% are replaced by column g
% if the new matrix is definite (signalled by p=[]);
% otherwise, the original L,d and
% a direction p of null or negative curvature are returned
%
% d contains diag(D) and is assumed positive
% Note that g must have zeros in other unit rows!!!
%
function [L,d,p]=ldlup(L,d,j,g);

p=[];

test=0;
if test,
    disp('enter ldlup')
    A=L*diag(d)*L';A(:,j)=g;A(j,:)=g';
end;

n=size(d,1);
I=1:j-1;K=j+1:n;
if j==1,
    v=zeros(0,1);
    del=g(j);
    if del<=n*eps,
        p=[1;zeros(n-1,1)];
        if test,
            A,p
            Nenner=abs(p)'*abs(A)*abs(p);
            if Nenner==0, indef1=0 ,else indef1=(p'*A*p)/Nenner, end;
            disp('leave ldlup at 1')
        end;
        return;
    end;
    w=g(K)/del;
    L(j,I)=v';
    d(j)=del;
    if test,
        A1=L*diag(d)*L',A
        quot=norm(A1-A,1)/norm(A,1),
        disp('leave ldlup at 3')
    end;
    return;
end;

% now j>1, K nonempty
LII=L(I,I);
u=LII\g(I);
v=u./d(I);
del=g(j)-u'*v;
if del<=n*eps,
    p=[LII'\v;-1;zeros(n-j,1)];
    if test,
        A,p
        indef1=(p'*A*p)/(abs(p)'*abs(A)*abs(p))
        disp('leave ldlup at 2')
    end;
    return;
end;
LKI=L(K,I);
w=(g(K)-LKI*u)/del;
[LKK,d(K),q]=ldlrk1(L(K,K),d(K),-del,w);
if isempty(q),
    % work around expensive sparse L(K,K)=LKK
    L=[L(I,:);
        v', 1,L(j,K);
        LKI,w,LKK];
    d(j)=del;
    if test,
        A1=L*diag(d)*L',A
        quot=norm(A1-A,1)/norm(A,1),
        disp('leave ldlup at 4')
    end;
else
    % work around expensive sparse L(K,K)=LKK
    L=[L(1:j,:);
        LKI,L(K,j),LKK];
    pi=w'*q;
    p=[LII'\(pi*v-LKI'*q);-pi;q];
    if test,
        indef2=(p'*A*p)/(abs(p)'*abs(A)*abs(p)),
        disp('leave ldlup at 5')
    end;
end;

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ldlrk1.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [L,d,p]=ldlrk1(L,d,alp,u)
% computes LDL^T factorization for LDL^T+alp*uu^T
% if alp>=0 or if the new factorization is definite
% (both signalled by p=[]);
% otherwise, the original L,d and
% a direction p of null or negative curvature are returned
%
% d contains diag(D) and is assumed positive
%
% does not work for dimension 0
%
function [L,d,p]=ldlrk1(L,d,alp,u);

test=0; % only for testing the routine
if test,
    disp('enter ldlrk1')
    A=L*diag(d)*L'+(alp*u)*u';
end;

p=[];
if alp==0, return; end;

n=size(u,1);
neps=n*eps;

% save old factorization
L0=L;d0=d;

% update
for k=find(u~=0)',
    del=d(k)+alp*u(k)^2;
    if alp<0 & del<=neps,
        % update not definite
        p=zeros(n,1);p(k)=1;
        p(1:k)=L(1:k,1:k)'\p(1:k);
        % restore original factorization
        L=L0;d=d0;
        if test,
            indef=(p'*(A*p))/(abs(p)'*(abs(A)*abs(p)))
            disp('leave ldlrk1 at 1')
        end;
        return;
    end;
    q=d(k)/del;
    d(k)=del;
    % in C, the following 3 lines would be done in a single loop
    ind=k+1:n;
    c=L(ind,k)*u(k);
    L(ind,k)=L(ind,k)*q+(alp*u(k)/del)*u(ind,1);
    u(ind,1)=u(ind,1)-c;
    alp=alp*q;
    if alp==0, break; end;
end;
if test,
    A1=L*diag(d)*L',A
    quot=norm(A1-A,1)/norm(A,1)
    disp('leave ldlrk1 at 2')
end;

end


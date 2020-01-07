% *************************************************************************
% Quantitative Macroeconomics: Final Project
% Marta Oliva Riera
% *************************************************************************
% COMPLEX VARIANT OF THE KRUSELL SMITH ALGORITHM
% Modifications on my code B2GE.m previously submitted for Project 3, now
% including aggregate shocks and trying to implement the Krusell-Smith algorithm
% Most changes are in the household problem, as I ran out of time
clear,clc


opt_det=false;          % 1=deterministic model
opt_nosr=false;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,
                        % 2=Markov chain with ny=2 (Krüger-Ludwig calibration)

% *************************************************************************
% SOLUTION FOR THE GENERAL EQUILIBRIUM MODEL: 
% *************************************************************************
global maxit tol df r nj ny replrate L R delta gridx tau kgrid nk alpha z

% Calibration:
func_calibr(opt_det,opt_nosr,opt_ny);

% made up guesses for psi
psi0 = [2;3];
psi1 = [0.5;1];

% approximation for the log of capital
lnk1=zeros(2,nk);
for j=1:nk
    for i=1:2
        lnk1(i,j)=psi0(i)+psi1(i)*log(kgrid(j));
    end
end
k1=exp(lnk1);
tau = func_pens(L,R,replrate);

% solution of household model
[gridx,gridsav,gridass,cfun,vfun] = func_hh(k1);


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function func_calibr(opt_det,opt_nosr,opt_ny)

global betta alpha tetta delta df r nj jr nx nk ny pi gridy pens sr epsi curv pini frac pop totpop grdfac maxit tol L  R replrate z piz kgrid

close all

r = 0.02;
rho = 0.04;
betta = 1/(1+rho);
tetta = 2;
delta = 0.05;   % depreciation rate of capital
alpha = 0.33;

nj=80;
jr=45;
replrate = 0.0;     % pension system replacement rate

nx=30;          % # of grid-points
curv=3.0;       % curvature of grid
grdfac=60;      % scaling factor of saving grid

tol = 1e-4;     % tolerance level
maxit = 100;    % maximum number of iterations on r
df = 0.1;       % "dampening factor" to update the guess on r

% deterministic income component:
% now wages and pensions will depend on the firm's problem and on the taxes
epsi=ones(nj,1);    % indicator for workers
if (jr<nj)
    epsi(jr+1:nj)=0.0;
end

% survival rates
if opt_nosr
    sr = ones(nj,1);
else
    mr = readfile([],'MR.txt',3);
    sr = 1.0-mr(21:21+nj-1,1);
end

% population and fraction living in year...
pop=zeros(nj,1);
pop(1)=100;
for jc=2:nj,
    pop(jc)=pop(jc-1)*sr(jc-1);
end
totpop=sum(pop);

% normalize population to one:
pop=pop/totpop;
totpop=1.0;
frac=pop./totpop;
L = sum(pop(1:jr));     % worker share of the population
R = sum(pop(jr+1:nj));  % retired share of the population

% # of income states
if (opt_det==1)        % deterministic model
    ny = 1;
    pini = 1.0;
    gridy = 1.0;
    pi = 1.0;
else                    % stochastic model
    
    if (opt_ny==1)      % markov chain with 5 states of the income shock
        % number of income shocks
        ny = 5;
        % transition probability
        rhoeta=0.98;
        % variance of "permanent" shock
        % taken from Campbell, Viceira, Ch. 7
        vareta=0.01;
        % Markov chain:
        [pi,gridy] = markovappr(rhoeta,sqrt(vareta),2,ny);
        
        % compute invariant distribution
        pini = 1/ny*ones(ny,1);
        for tc=1:100
            pini = pi'*pini;
        end
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else    % markov chain with two states - this is the one we want!
        
        % Alternative -- taken from Krüger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        % transition probability and variance
        rhoeta=0.97;
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % shock
        epsil=sqrt(vary/(4.0*rhoeta*(1.0-rhoeta)));
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
    end
    
end

% Aggregate technology shock with two states:
z = [1.03 0.97]';
piz = [0.95 0.05; 0.05 0.95];   % transition matrix

% Create a grid on aggregate capital: 
Keq = 8;    % equilibrium from project 3
nk = 5;
kgrid = linspace(0.5*Keq,1.5*Keq,nk);

end     % end function func_calibr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [gridx,gridsav,gridass,cfun,vfun] = func_hh(k1)
% solution to the household problem:
global alpha betta tetta delta tau nj nx nk ny pi gridy sr epsi curv grdfac replrate kgrid piz z

% grids and decisions rules:
gridx = zeros(nj,ny,nk,2,nx);
gridsav = zeros(nx,1);
gridass = zeros(nj,ny,nk,2,nx);
cfun = zeros(nj,ny,nk,2,nx);
vfun = zeros(nj,ny,nk,2,nx);
vpfun = zeros(nx,1);
vptrans = zeros(nj,ny,nk,2,nx);

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% income states
for yc=1:ny
    for kc=1:nk
        for zc=1:2
            % current wages and interest rate:
            wage=(1-alpha)*kgrid(kc)^(alpha)*z(zc);
            ret=alpha*kgrid(kc)^(alpha-1)*z(zc)-delta;
            
            % cash-on-hand grid at nj:
            inc = epsi(nj)*wage*(1-tau)*gridy(yc)+(1-epsi(nj))*replrate*wage*(1-tau);

            % in case of no pension system, assume some minimum cash on hand:
            minx=max(inc,sqrt(eps));
            maxx=gridsav(nx)*(1.0+ret)+inc;
            gridx(nj,yc,kc,zc,:)=linspace(minx,maxx,nx);

            % Final period Consumption function, asset holdings, value function, including derivative
            cfun(nj,yc,kc,zc,:)=gridx(nj,yc,kc,zc,:);
            gridass(nj,yc,kc,zc,:)=(gridx(nj,yc,kc,zc,:)-inc)/(1+ret);
            vfun(nj,yc,kc,zc,:)=U(cfun(nj,yc,kc,zc,:));
            vpfun(:)=MUc(cfun(nj,yc,kc,zc,:));
            vptrans(nj,yc,kc,zc,:)=vpfun.^(-1.0/tetta);
        end
    end
end

% Iterate Backwards
for jc=nj-1:-1:1
    for yc=1:ny
        for kc=1:nk
            for zc=1:2
                K11=k1(zc,kc);
                
                %current wage and interest rate
                wage=(1-alpha)*kgrid(kc)^(alpha)*z(zc);
                ret=alpha*kgrid(kc)^(alpha-1)*z(zc)-delta;
                
                % interpolate K1 on the grid of aggregate capital
                A = repmat(K11,[1 length(kgrid)]);
                [xx, stk]=min(abs(A-kgrid));
                K1=kgrid(stk);
                
                for xc=2:nx,
                    vp=zeros(2,2);                
                    for ycc=1:ny
                        for zcc=1:2
                            % income tomorrow:
                            w1=(1-alpha)*K1^alpha*z(zcc);
                            ret1=alpha*K1^(alpha-1)*z(zcc)-delta;
                            incp1=epsi(jc+1)*w1*(1-tau)*gridy(ycc)+(1-epsi(jc+1))*replrate*w1*(1-tau);

                            % Maximum cash on hand tomorrow:
                            % in case of zero savings and no pension system assume some
                            % minimum cash on hand
                            cah=max(sqrt(eps),incp1+(1.0+ret1)*gridsav(xc));

                            % Interpolate derivative of value function
                            if ( cah<gridx(jc+1,ycc,stk,zcc,1))
                                disp('how can this be?')
                            end
                            if ( cah>gridx(jc+1,ycc,stk,zcc,nx) )
                                % if out of bounds simply set it to decision at nx:
                                vptr = vptrans(jc+1,ycc,stk,zcc,nx);
                            else
                                vptr = interp1(squeeze(gridx(jc+1,ycc,stk,zcc,:)),squeeze(vptrans(jc+1,ycc,stk,zcc,:)),cah);
                            end
                            vp(ycc,zcc)=vptr.^(-tetta); %value function
                        end %for zcc
                    end %for ycc
                    
                PI = pi(yc,:)'*piz(zc,:);
                V = vp.*PI;

                % Euler equation: RHS
                expvp=betta*sr(jc)*(1.0+ret)*sum(sum(V));

                % consumption
                cfun(jc,yc,kc,zc,xc)=invut(expvp);

                % endogenous x-grid:
                gridx(jc,yc,kc,zc,xc)=gridsav(xc)+cfun(jc,yc,kc,zc,xc);
                end %for xc
                    
            % income (wages and pensions) in current period/age:
            inc=epsi(jc)*wage*(1-tau)*gridy(yc)+(1-epsi(jc))*replrate*wage*(1-tau);

            % decision at minx
            % notice: correction required for welfare calculation
            % the above is actually slightly inefficient because xmin
            % can be explicitly computed, then gridsav would be age and
            % state dependent.
            minx=max(inc,sqrt(eps));
            if (minx<gridx(jc,yc,kc,zc,2))
                gridx(jc,yc,kc,zc,1)=minx;
            else    % set it to some arbitrary fracion of x(2)
                gridx(jc,yc,kc,zc,1)=0.9*gridx(jc,yc,kc,zc,2);
            end

            % Compute optimal consumption and leisure for minx
            cfun(jc,yc,kc,zc,1)=gridx(jc,yc,kc,zc,1);

            % assets at all xc:
            gridass(jc,yc,kc,zc,:)=(gridx(jc,yc,kc,zc,:)-inc)/(1+ret);

            % Update vfun and vpfun
            vpfun(:)=MUc(cfun(jc,yc,kc,zc,:));
            vptrans(jc,yc,kc,zc,:)=vpfun(:).^(-1.0/tetta);

            % Calculate value function
            for xc=1:nx

                v=zeros(2,2);
                for ycc=1:ny
                    for zcc=1:2
                        % income tomorrow:
                        w1=(1-alpha)*K1^alpha*z(zcc);
                        ret1=alpha*K1^(alpha-1)*z(zcc)-delta;
                        incp1=epsi(jc+1)*w1*gridy(ycc)*(1-tau)+(1-epsi(jc+1))*w1*replrate*(1-tau);

                        % cah tomorrow
                        cah=max(sqrt(eps),incp1+(1.0+ret1)*gridsav(xc));

                        % this should never be the case:
                        if ((cah+0.0001)<gridx(jc+1,ycc,stk,zcc,1))
                            warning('How can this be ?');
                        end
                        % linear interpolation:
                        v(ycc,zcc)=func_intp(squeeze(gridx(jc+1,ycc,stk,zcc,:)),squeeze(vfun(jc+1,ycc,stk,zcc,:)),cah);
                    end % for zcc
                end    % end for ycc

                PI=pi(yc,:)'*piz(zc,:);
                V=v.*PI;
 
                % update value function
                expv=sum(sum(V));
                vfun(jc,yc,kc,zc,xc)=U(cfun(jc,yc,kc,zc,xc))+betta*sr(jc)*expv;
                
                end    % end for xc
            end %for zc
        end % for kc
    end %for yc
end    % end for jc


% ---------------------------------------------------------------------
    function fv = func_intp(x,func,xp)
        
        
        n = length(x);
        if ( xp>x(n) )
            % fv = func(n);
            fv=func_extrapol(x(n-1),x(n),func(n-1),func(n),xp);
        elseif (xp<x(1))
            % fv = func(1);
            fv=func_extrapol(x(1),x(2),func(1),func(2),xp);
        else
            fv = interp1(x,func,xp);
        end
        
    end
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
    function y=func_extrapol(x1,x2,y1,y2,x)
        
        % simple linear extrapolation
        
        m = (y2-y1)/(x2-x1);
        y = y1 + m*(x-x1);
        
    end
% ---------------------------------------------------------------------

end     % end function func_hh
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [Phi,PhiAss,ass]=func_aggr(gridx,gridsav,cfun,gridass,wage,ret,tau)
% aggregation and cross-sectional measure
global nj nx ny pi gridy pens sr epsi pini frac totpop replrate

% Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx);          % distribution of assets conditional by age and shock
PhiAss = zeros(nx,1);             % distribution of assets

% Distribution of newborns over cash at hand
for yc=1:ny
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*wage*(1-tau)*gridy(yc)+(1-epsi(1))*replrate*wage*(1-tau);
    
    % initial cash-on-hand:
    cahini=inc;
    
    [vals,inds]=basefun(gridx(1,yc,:),cahini,nx);
    Phi(1,yc,inds(1))=vals(1)*pini(yc)*frac(1);
    Phi(1,yc,inds(2))=vals(2)*pini(yc)*frac(1);
end

for jc=2:nj
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
        for yc=1:ny
            for ycc=1:ny
                
                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*wage*(1-tau)*gridy(ycc)+(1-epsi(jc))*replrate*wage*(1-tau);
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+ret)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ycc,:),cah,nx);
                
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
            end    
        end    
    end    
   
    for xc=1:nx
        for yc=1:ny
            for xcc=1:nx
                for ycc=1:ny
                    % transfer distribution:
                    Phi(jc,ycc,xcc)=Phi(jc,ycc,xcc)+Phi(jc-1,yc,xc)*TT(ycc,xcc,yc,xc)*sr(jc-1);
                end
            end
        end
    end
    
end    % end for jc

% Check that for each country distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) )
    beep; beep; beep;
    warning('distribution fucked up');
end

% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 )
    beep; beep; beep;
    warning('grid too small -- increase your grid');
    pause
end

ass=0.0;
cons=0.0;
lab=0.0;
re=0.0;

% aggregation
for jc=1:nj
    for yc=1:ny
        for xc=1:nx
            PhiAss(xc)=PhiAss(xc)+Phi(jc,yc,xc);
            
            % asset holdings = capital stock in general equilibrium
            ass=ass+totpop*Phi(jc,yc,xc)*gridsav(xc);
            
            cons=cons+totpop*Phi(jc,yc,xc)*cfun(jc,yc,xc);
            
            lab=lab+totpop*Phi(jc,yc,xc)*gridy(yc)*epsi(jc);
            re=re+totpop*Phi(jc,yc,xc)*gridy(yc)*(1.0-epsi(jc));
        end
    end
end


% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookup(grid_x,x,0);
        
        if ( (i+1)>nx)
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0)
            inds(1)=1;
            inds(2)=1;
            vals(1)=1.0;
            vals(2)=0.0;
        else
            inds(1)=i;
            inds(2)=i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2)=( x-grid_x(i) )/dist;
            vals(1)=( grid_x(i+1)-x )/dist;
        end
        
    end 	% end function basefun
% ---------------------------------------------------------------------




end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [pini,pi,gridy]=mchain(rhoeta,epsil)

% Transition Probabilities
pi=rhoeta*ones(2,2);
pi(1,2)=1.0-rhoeta;
pi(2,1)=1.0-rhoeta;

% Initial Distribution
pini=0.5*ones(2,1);

gridy=zeros(2,1);
gridy(1)=exp(1.0-epsil);
gridy(2)=exp(1.0+epsil);
gridy=2.0*gridy/(sum(gridy));

end  % end function mchain
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function u = U(c)
global tetta

% utility
if (abs(tetta-1-0)<sqrt(eps))
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end
end     % end function U
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function muc=MUc(c)
global tetta

% maringal utility
muc = c.^(-tetta);
end     % end function MUc
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function invut=invut(marg)
global tetta

% invert utility for c
invut=marg.^(-1.0/tetta);
end     % end function invut
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grd = makegrid(x1,x2,n,c)
% makes curved grid according to curvature parameter c
scale=x2-x1;
grd(1)=x1;
grd(n)=x2;
for i=2:n-1
    grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
end
end     % end function makegrid
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);     % tax on income

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [mpk,Y] = func_mpk(ass, L)

global alpha

Y = ass.^alpha * L.^(1-alpha);
ky = ass./Y;
mpk = alpha * ky.^(-1);

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



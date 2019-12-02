% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Exogenous grid solution method for the household problem
% Marta Oliva Riera
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear, clc;
global nj ny 

tic

opt_det=false;          % 1=deterministic model
opt_nosr=false;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,
                        % 2=Markov chain with ny=2 (Krüger-Ludwig calibration)

% -------------------------------------------------------------------------
% SOLUTION
func_calibr(opt_det,opt_nosr,opt_ny);
[gridx,gridsav,gridass,cfun,vfun] = EXOGM_hh;

toc

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function func_calibr(opt_det,opt_nosr,opt_ny)

global xmax xmin betta tetta r nj jr nx ny pi gridy netw pens sr epsi curv pini frac pop totpop grdfac

close all

r = 0.04;
rho = 0.04;
betta = 1/(1+rho);
tetta = 2;

nj=80;
jr=45;

nx=30;              % # of grid-points
curv=3.0;           % curvature of grid
xmax = 30;          % scaling factor of saving grid
xmin = sqrt(eps);

% deterministic income component:
netw=1.0;
pens=0.4;
epsi=ones(nj,1);
if (jr<nj),
    epsi(jr+1:nj)=0.0;
end;

% survival rates
if opt_nosr,
    sr = ones(nj,1);
else
    mr = readfile([],'MR.txt',3);
    sr = 1.0-mr(21:21+nj-1,1);
end;

% population and fraction living in year...
pop=zeros(nj,1);
pop(1)=100;
for jc=2:nj,
    pop(jc)=pop(jc-1)*sr(jc-1);
end;
totpop=sum(pop);

% normalize population to one:
pop=pop/totpop;
totpop=1.0;
frac=pop./totpop;

% # of income states
if (opt_det==1),        % deterministic model
    ny = 1;
    pini = 1.0;
    gridy = 1.0;
    pi = 1.0;
else                    % stochastic model
    
    if (opt_ny==1)      % markov chain with 5 states of the income shock - not our case for now
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
        for tc=1:100,
            pini = pi'*pini;
        end;
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else    % markov chain with two states
        
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
    end;
    
end;


end     % end function func_calibr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [gridx,gridsav,gridass,cfun,vfun] = EXOGM_hh

global vpfun betta tetta r nj nx ny pi gridy netw pens sr epsi curv xmax xmin

disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nj,ny,nx);
gridass = zeros(nj,ny,nx);
cfun = zeros(nj,ny,nx);
vfun = zeros(nj,ny,nx);
vpfun = zeros(nx,1);

% cash on hand grid:
gridx=makegrid(xmin,xmax,nx,curv);
gridx=gridx';

% for each income state:
for yc=1:ny 
    % Final period Consumption function, asset holdings, value function, including derivative
    cfun(nj,yc,:)=gridx;
    vfun(nj,yc,:)=U(cfun(nj,yc,:));
end;

% Iterate Backwards: 
for jc=nj-1:-1:1,           % on age (j)
    for yc=1:ny             % for every income shock state (ny)
        vpfun(:) = MUc(cfun(jc+1,yc,:));
        for xc=1:nx,        % for every point in the grid of x (nx)
            % check if the borrowing constraint binds:
             mincons=gridx(xc);
             mu = foc3(mincons,gridx(xc),jc,yc);
             if (mu>=0.0),
                 cfun(jc,yc,xc)=mincons;
             else
                 [cfun(jc,yc,xc),fval] = fzero('foc3',1,[],gridx(xc),jc,yc);
                 cfun(jc,yc,xc)=max(cfun(jc,yc,xc),0.00001);
             end;
             
            if cfun(jc,:,xc) > gridx(jc,:,xc) + phi/(1+r);
                cfun(jc,:,xc) = gridx(jc,:,xc) + phi/(1+r);
            end
            
            end;
         inc=epsi(jc)*netw*gridy(yc)+(1-epsi(jc))*pens;
         gridass(jc,yc,:)=(gridass(jc+1,yc,:)-inc+cfun(jc,yc,xc))/(1+r);
     end
end
            
end     % end function EXOGM_hh
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function fv = func_intp3(x,func,xp)

        n = length(x);
        if ( xp>x(n) ),
            % fv = func(n);
            fv=func_extrapol3(x(n-1),x(n),func(n-1),func(n),xp);
        elseif (xp<x(1)),
            % fv = func(1);
            fv=func_extrapol3(x(1),x(2),func(1),func(2),xp);
        else
            fv = interp1(x,func,xp);
        end;
        
  end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function y=func_extrapol3(x1,x2,y1,y2,x)
        
        % simple linear extrapolation
        
        m = (y2-y1)/(x2-x1);
        y = y1 + m*(x-x1);
        
end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function fval = foc3(cons,x,t,yc)

global betta tetta r g gridx vpfun epsi probepsi ne nx sr

vpp1 = evalvp3(cons,x,t,yc);
margu = MUc(cons);
fval = margu - betta*sr(t)*(1.0+r) * vpp1;

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function vpp1 = evalvp3(cons,x,t,yc)

global betta tetta r g gridx vpfun epsi netw gridy pens pi ny nx 
vpp1 = zeros(ny,1);
for ec=1:ny,
    xp1 = (x-cons)*(1+r)+epsi(t)*netw*gridy(ec)+(1-epsi(t))*pens;
    vpp1(ec) = (func_intp3(gridx,vpfun,xp1));
end;
vpp1 = sum(vpp1'.*pi(yc,:));

end
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
if (abs(tetta-1-0)<sqrt(eps)),
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end;
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

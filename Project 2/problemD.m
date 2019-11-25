% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Quantitative Macroeconomics: Project 2
% Marta Oliva Riera
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clear, clc
global betta tetta r g gridx vpfun epsi probepsi ne nx min_cons

% PROBLEM D: FINITE TIME HORIZON
% -------------------------------------------------------------------------
% SETTINGS
tol=1e-4;
min_cons=1.0e-08;

nj = 80;    % maximum age and last period, like in Project 1

% parameters
r = 0.02;
rho = 0.03;
g = 0.01;
tetta = 1;
betta = 1/(1+rho);

% grid
nx=50;              % # of grid-points
curv=3.0;           % curvature of grid
xmax = 30;          % scaling factor of saving grid
xmin = sqrt(eps);
gridx=makegrid(xmin,xmax,nx,curv);
gridx=gridx';

% income shocks
ne = 7;
varepsi = 0.01;
muepsi = -varepsi/2;
[epsi,probepsi] = qnwnorm(ne,muepsi,varepsi);
mat=[[1:ne]',epsi,probepsi];
epsi=exp(epsi);

if (abs(sum(epsi.*probepsi)-1.0)>sqrt(eps)),
    error('random numbers fucked up');
end;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SOLUTION 
% initial guess
cons = gridx;                 % consumption
% vpfun = margutil(cons);     % derivative of value function
xT1 = 0;    % future cash on hand in the last period
minconsT = gridx(end);
muT = foc(minconsT, xT1);
if (muT>=0.0),           % if mu >= 0 then the constraint is binding, so x=c
    cons(nj)=minconsT;
else,                   % if mu < 0 we have an interior solution: interpolate on the FOC
    [cons(nj),fval] = fzero('foc',cons(nj),[],gridx(nj));   
    cons(nj)=max(cons(nj),min_cons);
end;

% Iteration
for it = nj-1:-1:1
    % check binding constraint:
    mincons=gridx(it);
    mu = foc(mincons,gridx(it));
    if (mu>=0.0),           % if mu >= 0 then the constraint is binding, so x=c
        cons(it)=mincons;
    else,                   % if mu < 0 we have an interior solution: interpolate on the FOC
        [cons(it),fval] = fzero('foc',cons(it),[],gridx(xc));   
        cons(it)=max(cons(it),min_cons);
    end;
    % update vpfun
    vpfun = margutil(cons);
end;    

% -------------------------------------------------------------------------

figure;
plot(gridx,cons,'LineWidth',2);
xlabel('x');
ylabel('c');
title('consumption policy');

figure;
plot(gridx(3:end),vpfun(3:end),'LineWidth',2);
xlabel('x');
ylabel('V(x)');
title('Value function');

% -------------------------------------------------------------------------
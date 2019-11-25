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
cons = [];
cons(:, nj) = gridx;       % consumption(T) = cash at hand(T)

for t=nj-1:-1:1,
    vpfun = margutil(cons(:,t+1));
    for xc=1:nx,    % for each point in the grid
        % check binding constraint:
        mincons=gridx(xc);
        mu = foc(mincons,gridx(xc));
        if (mu>=0.0),           % if mu >= 0 then the constraint is binding, so x=c
            cons(xc, t)=mincons;
        else,                   % interior solution
            [cons(xc, t),fval] = fzero('foc',cons(xc, t),[],gridx(xc));   
            cons(xc, t)=max(cons(xc, t),min_cons);
        end;
    end;
end

% -------------------------------------------------------------------------

figure;
plot(gridx,cons(:,nj),'LineWidth',2);
xlabel('x');
ylabel('c');
title('consumption policy (t=T)');


% -------------------------------------------------------------------------
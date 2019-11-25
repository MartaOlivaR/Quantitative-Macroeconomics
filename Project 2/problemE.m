% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% project 2
% solution of household problem for T = \infty

clear all
close all

global betta tetta r g gridx vpfun epsi probepsi ne nx min_cons

% -------------------------------------------------------------------------
% SETTINGS
maxit = 100; 
tol=1e-4;
nt=1100;    % periods for simulation
dt=100;     % periods discarded
min_cons=1.0e-08;

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
% [epsi,probepsi] = qnwnorm(ne,muepsi,varepsi);
% mat=[[1:ne]',epsi,probepsi];
% save rn.txt mat -ascii -double -tabs
mat=load('rn.txt');
epsi=mat(:,2);
probepsi=mat(:,3);
epsi=exp(epsi);
if (abs(sum(epsi.*probepsi)-1.0)>sqrt(eps)),
    error('random numbers fucked up');
end;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SOLUTION 
% initial guess
cons = gridx;               % consumption
vpfun = margutil(cons);     % derivative of value function

% iteration
for it=1:maxit,
    disp(['iteration # ', num2str(it)]);
    consm1=cons;
    
    for xc=1:nx,    % for each point in the grid
        % check binding constraint:
        mincons=gridx(xc);
        mu = focE(mincons,gridx(xc));
        if (mu>=0.0),           % if mu >= 0 then the constraint is binding, so x=c
            cons(xc)=mincons;
        else,   
            [cons(xc),fval] = fzero('focE',cons(xc),[],gridx(xc));   
            cons(xc)=max(cons(xc),min_cons);
        end;
%          if (cons(xc)>gridx(xc))
%              cons(xc)=gridx(xc);
%          end;
    end;    

    % update vpfun
    vpfun = margutil(cons);
    
    % check convergence
    dist=cons-consm1;
    maxdist=max(abs(dist));
    if (maxdist<tol),
        disp(['I am so happy! The thing has converged in iteration ', num2str(it)]);
        break;
    else,
        disp(['current maximum norm is ', num2str(maxdist)]);
    end;
end;
if (it==maxit),
    warning('increase # of iters for specified tolerance');
end;

% -------------------------------------------------------------------------

figure;
plot(gridx,cons,'LineWidth',2);
xlabel('x');
ylabel('c');
title('consumption policy');

%[gridx,cons]



%% Plots to compare the value function and consumption
plot(vpfun(5:end))
plot(cons)

%%
function fval = focE(cons,x)

global betta tetta r g gridx vpfun epsi probepsi ne nx 

vpp1 = evalvpE(cons,x);
margu = margutil(cons);
fval = margu - (betta * (1+g)^(-tetta) * (1+r) * vpp1);

end
% -------------------------------------------------------------------------
function vpp1 = evalvpE(cons,x)
% For each x not on the grid, interpolate its corresponding inverse value
% function (A). Then invert A to get the Value function and take expectations 
% over the possible states of the income shock
global betta tetta r g gridx vpfun epsi probepsi ne nx 

vpp1 = zeros(ne,1);
for ec=1:ne,
    xp1 = (x-cons)*(1+r)/(1+g)+epsi(ec);
    App1(ec) = func_intp(gridx,cons,xp1);
end;
vpp1 = inv(App1)
vpp1 = sum(vpp1.*probepsi); 

end
% -------------------------------------------------------------------------
 function fv = func_intp(x,func,xp)
 chkout=false;
    n = length(x);
    if ( xp>x(n) ),
        % fv = func(n);
        fv=func_extrapol(x(n-1),x(n),func(n-1),func(n),xp);
    elseif (xp<x(1)),
        % fv = func(1);
        fv=func_extrapol(x(1),x(2),func(1),func(2),xp);
    else
        fv = interp1(x,func,xp);
    end;
    
 end
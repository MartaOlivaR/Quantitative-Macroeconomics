% *************************************************************************
% Quantitative Macroeconomics: Final Project
% Marta Oliva Riera
% *************************************************************************
% SIMPLE VARIANT OF THE KRUSELL SMITH ALGORITHM
%% (2) Simulate a first order difference equation in logs:
%--------------------------------------------------------------------------
clear, clc
rng(42);

% Define the parameter values:
g = 0;              % growth rate of technology level
lambda = 0.5;       % relative length of working in j=2 / labour productivity in j=2
alpha = 0.3;        % capital elasticity in Cobb Douglas production function
beta = 0.99^40;     % depreciation rate
tau = 0;          % tax level
T = 50000;          % number of periods

% Set up the iid shocks:
% aggregate shocks, with two possible states each (mean +- std):
lnzeta = [-sqrt(40*0.02^2); sqrt(40*0.02^2)];     % shock to wages
zeta = exp(lnzeta);
lnrho = [-sqrt(40*0.08^2);sqrt(40*0.08^2)];       % shock to rate of return
rho = exp(lnrho);

% idiosyncratic shock to wages, discretized with 11 nodes:
ne = 11;                                   % grid size
sig2e = 40*0.15^2;                         % variance log(eta)
[lneta,probeta] = qnwnorm(ne,0,sig2e);     % gaussian quadrature
eta = exp(lneta);

% Compute phi, as an expectation given the different possible states of the shocks:
Phi = zeros(2,11);

for i = 1:2
    for j = 1:11
        Phi(i,j) = 1/(1+((1-alpha)/(alpha*(1+lambda)*rho(i)))*(lambda*eta(j)+tau*(1+lambda*(1-eta(j)))));
    end
end

phi = sum(Phi(1,:)*probeta*0.5)+sum(Phi(2,:)*probeta*0.5);  
s = (beta*phi)/(1+beta*phi);    % savings

% Steady state of capital:
lnk0 = (log(s) + log(1-tau) - log(1+lambda) + log(1-alpha) )/(1-alpha);
k0 = exp(lnk0);

% Simulation of the first order difference equation for the log of capital:
lnk = zeros(T,1);
state = zeros(T,1);
lnz = zeros(T,1);
lnk(1) = lnk0;    % steady state as the initial value   

for t = 1:T
   state(t) = randi(length(zeta));     % randomly select boom or recession in each t
   lnz(t) = lnzeta(state(t));          % realizations of the zeta shock for every t
   z(t) = exp(lnz(t));
   lnk(t+1) = log(s) + log(1-tau) + log(1-alpha) - log(1+lambda) + lnz(t) + alpha*lnk(t);
end

k = exp(lnk);
%plot(k)

%% (3) Krusell Smith algorithm:
%--------------------------------------------------------------------------
% Some more parameters which need to be defined:
upsilon = 1;        % technology level
epsilon = 0.0001;   % tolerance level for convergence
omega = 0.2;        % coefficient to update the guess of psi in the iterations
maxit = 100;        % maximum iterations in the algorithm
n = 5;              % number of grid points for capital

% Create a grid for capital, with equally spaced points:
kgrid = linspace(log(0.5)+lnk0,log(1.5)+lnk0,n);

% (a) Theoretical values of psi:
psi1 = [alpha;alpha];
psi0 = log(s) + log(1-tau) + log(1-alpha) + lnzeta - log(1+lambda);


% (b) Simple Krusell-Smith algorithm
for m = 1:maxit
    
    k = [kgrid;kgrid];              % capital grid for both aggregate states
    lnk1 = psi0*ones(1,5)+psi1.*k;  % approximated next period ln(k) using psi
    k1 = exp(lnk1);
    psi = [psi0 psi1];
    
    for i = 1:2
        % current period wages
        w(i,:) = (1-alpha)*upsilon*exp(k(i,:)).^alpha*zeta(i);
        % next period interest rate and wages
        R1(i,:) = alpha*k1(i,:).^(alpha-1)*zeta(i)*rho(i);
        w1(i,:) = (1-alpha)*upsilon*k1(i,:).^alpha*zeta(i); 
    end
    
    
    % Household problem:
    % expectation terms in the assets variable
    b = (w1*R1'.^(-1))*0.5*0.5;
    B = sum(b);
    c = eta.*probeta;
    C = sum(c');
    % assets:
    a = (beta*(1-tau)*w - B'*(lambda*(1-tau)*C - tau*(1+lambda)))/(1+beta);
    for j=1:n       % adding a non-negativity constraint for assets (otherwise you get complex numbers)
        for i = 1:2
            if a(i,j) < 0 
                a(i,j) = 0.00001;
            end
        end
    end
    
    s = a./((1-tau)*w);
    c1 = (1-tau)*w-a;
    %c2 = a*R1' + lambda*eta*w1*(1-tau)+tau*w1*(1+lambda)
    
    
    % Simulation:
    lnkS = zeros(1,T);
    ss = 3;
    lnkS(1) = kgrid(ss);     % starting in steady state
    for t = 2:T
        lnk1S = log(s(state(t),ss)) + log(1-tau) + log(1-alpha) - log(1+lambda) + lnz(t) + alpha*lnkS(t-1);
        lnkS(t) = lnk1S;
    end
    kS = exp(lnkS);
        
    % Regression:
    bcount = 1; rcount = 1;
    kB = []; k1B = [];
    kR = []; k1R = [];
    
    for t = 501:T-1     % getting rid of the first 500 observations
        if lnz(t) >= 0  % if in a boom, save the simulated lnk in kB
            kB(bcount) = lnkS(t);
            k1B(bcount) = lnkS(t+1);
            bcount = bcount+1;
        else            % if in a recession, save the simulated lnk in kR
            kR(rcount) = lnkS(t);
            k1R(rcount) = lnkS(t+1);
            rcount = rcount+1;
        end
    end
    
    k1B = k1B';
    k0B = [ones(length(kB),1) kB'];
    k1R = k1R';
    k0R = [ones(length(kR),1) kR'];
    
    % regressions for each aggregate state (boom or recession):
    psiB = (k0B'*k0B)^(-1)*k0B'*k1B;
    psiR = (k0R'*k0R)^(-1)*k0R'*k1R;
    PSI = [psiR'; psiB'];   % saving the coefficients in the same way as the initial ones
    
    
    % Check convergence:
    if abs(psi - PSI) < epsilon     % if convergence is found, stop the loop
       display('Convergence reached in ' + string(m) + ' iterations!')
       break
    else        % if there's no convergence, update the guesses on psi and run the loop again
        psi0(1) = omega*psiR(1) + (1-omega)*psi0(1);
        psi0(2) = omega*psiB(1) + (1-omega)*psi0(2);
        psi1(1) = omega*psiR(2) + (1-omega)*psi1(1);
        psi1(2) = omega*psiB(2) + (1-omega)*psi1(2);     
    end
    
end

%% (d) Welfare comparison:
% some more parameters
theta = 1;
Beta = beta/(1+beta);   % relative importance of the second period of life

% simulation to obtain expected c2
for t = 1:T-1
    % wage in the current period and consumption of the young generation at t
    w = (1-alpha)*kS(t)^alpha*z(t);
    c1S = (1-tau)*w - kS(t+1);
    
    % interest rate and wages in the next period
    R1 = alpha*kS(t+1)^(alpha-1)*rho*zeta';
    w1 = (1-alpha)*upsilon*kS(t+1)^alpha*zeta;
    
    % expected consumption of the currently young (j=1) next period (when old, j=2)
    d = eta*w1';
    D = sum(d(:,1)*probeta(1)*0.5) + sum(d(:,2)*probeta(2)*0.5);
    E = sum(R1(1))*0.5 + sum(R1(2))*0.5;
    Ec2 = kS(t+1)*E + lambda*(1-tau)*D + tau*(1+lambda)*sum(w1)*0.5;
    
    U2(t) = log(Ec2);                           % expected utility when old
    EU(t) = (1-Beta)*log(c1S) + Beta*U2(t);     % expected lifetime utility for the generations born at every t
            
end

% Average expected utility over the whole simulation (excluding the first
% 500 observations):
avgEU = 1/(T-500)*sum(EU(501:end));

%% Consumption Equivalent Variation:
% Saving the results as a txt file to be able to compute CEV (run all of it
% twice with the different settings of tau for it to work):
if tau == 0
    save v0.txt avgEU -ascii -double
elseif tau == 0.1
    save v1.txt avgEU -ascii -double
end

V0 = importdata('v0.txt');     % Utility before the tax
V1 = importdata('v1.txt');     % Utility after the tax is introduced
g = exp((V1-V0)/beta)-1;        % CEV


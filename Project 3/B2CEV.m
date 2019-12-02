% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Consumption equivalent variation:
% Marta Oliva Riera
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

clear, clc

v0 = importdata('v0.txt');  % obtained from B2GE with different settings
v1 = importdata('vT.txt');  % obtained from B2GE(v1.txt) or B2CEVpartialeq(vT.txt)

theta = 2;

% CEV:
g = ((v1./v0).^(1/(1-theta)))-1

% Author: Matt Dvorsky

%% Clear Workspace
clc;
clear;
close all;

%% Inputs
er = {4 - 0.001j, 2 - 0.001j};
ur = {};
thk = {2, 0.1};

%% Create nLayer Object
NL = nLayerCircular(0, 5, waveguideBand="Ka_TE01", modeSymmetryAxial="TE");

%% Plot
figure;
nLayerViewer(er, ur, thk, NL, [32, 40]);


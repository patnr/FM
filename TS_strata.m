% Emulate layers in the ground and hence the wave front of an
% earthquake.

%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

m  = 100;
n  = 100;
dy = 1;
dx = 1;
Ly4 = m*dy;
Lx4 = n*dx;
Dxyz = [dx dy];

nStrata = 4;
F = ones(m,n) + 10 * repmat(ceil(linspace(eps,nStrata,m)'),1,n);

SPs = [30;40];


%% Calculation
T4 = fm(F,SPs,Dxyz,struct('implementation','C++','order',2));


%% Plot
figure(1);
% imshow(T4,[],'Init','fit','XD',[1 Lx4],'YD',[1 Ly4]);
contourf(flipud(T4)); box off;
title('Wave fronts in layered medium');
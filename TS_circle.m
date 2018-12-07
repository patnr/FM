% Circle (disc) test of fm2d with error measurements

%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

% Domain dimensions
Lx = 80; 
Ly = 100;

% Number of grid nodes in each direction
n = 20;
m = 200;

% An interesting experiment is to uncomment the 'showProgress' line in
% fm2d.m, which will result in a "video" of the fast marching process.
% Must use low m,n. Like 20.

% Dxyz
dx = Lx/n;
dy = Ly/m;
Dxyz = [dx dy];

% Speed map
F = ones(m,n);
% F = ones(m,n) + 4*rand(m,n);
% load('data/F.mat','F');

% Source points
SPs = [0 0; 62 93]';


%% Calculate Fast-marching solutions
T1 = fm(F,SPs,Dxyz,'implementation','Matlab','order',2);
T2 = fm(F,SPs,Dxyz,'implementation','C++','order',2);


%% Calculate exact answer
% Calculate grid node positions
[xx yy] = fmMeshGrid([n m], Dxyz);

% Calculate one time/distance map for each SP
Texa = zeros(m,n,size(SPs,2));
for iter=1:size(SPs,2)
	Texa(:,:,iter) = ...
		sqrt((xx-SPs(1,iter)).^2 + (yy-SPs(2,iter)).^2);
end
Texa = min(Texa,[],3);


%% Calculate error
rel_errors = (T1-Texa)./(Texa);
per_pixel_rel_error = sum(abs(rel_errors(:)))/(m*n);

display(per_pixel_rel_error);


%% Plots
% Get max values to enforce coherent colormaps
minT = min([T1(:)' Texa(:)']);
maxT = max([T1(:)' Texa(:)']);
minD = min(rel_errors(:));
maxD = max(rel_errors(:));


figure(1);

subplot(2,2,1);
imshow(T1,[minT maxT],...
	'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Fast marching solution');

subplot(2,2,2);
imshow(T2,[],'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Fast marching solution 2');

subplot(2,2,3);
imshow(Texa,[minT maxT],...
	'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Exact solution');

colorbar('Position', [.01 .1 .05 .8]);

subplot(2,2,4);
imshow(rel_errors,[minD maxD],...
	'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Relative errors');

colorbar('Position', [.87 .1 .05 .8]);

colormap(jet());
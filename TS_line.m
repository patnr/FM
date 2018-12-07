% One dimensional fm2d

%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

% Number of grid nodes in each direction
n = 100;
m = 1;

% Domain dimensions
Lx = 100; 
Ly = 20;

% Dxyz
dx = Lx/n;
dy = Ly/m;
Dxyz = [dx dy];

% Speed map
F = ones(m,n);
% F = ones(m,n) + 4*rand(m,n);

% Source points
SPs = [0 10]';


%% Calculate Fast-marching solutions
T1 = fm(F,SPs,Dxyz,'impl','C++');
T2 = fm(F,SPs,Dxyz,'impl','Matlab');
% T2 = fm(F,SPs,Dxyz,'impl','noHeap');


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
figure(1);

subplot(2,2,1);
imshow(T1,[],'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Fast marching solution');

subplot(2,2,2);
imshow(T2,[],'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Fast marching solution 2');

subplot(2,2,3);
imshow(Texa,[],'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Exact solution');

colorbar('Position', [.01 .1 .05 .8]);

subplot(2,2,4);
imshow(rel_errors,[],'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Relative errors');

colorbar('Position', [.9 .1 .05 .8]);

colormap(hot());
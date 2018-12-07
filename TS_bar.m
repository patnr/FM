% One dimensional fm3d

%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

% Number of grid nodes in each direction
n = 1;
m = 1;
o = 100;

% Domain dimensions
Lx = 20; 
Ly = 20;
Lz = 100;

% Dxyz
dx = Lx/n;
dy = Ly/m;
dz = Lz/o;
Dxyz = [dx dy dz];

% Speed map
F = ones(m,n,o);
% F = ones(m,n) + 4*rand(m,n);

% Source points
SPs = [10 10 0]';


%% Calculate Fast-marching solutions
T1 = fm(F,SPs,Dxyz,struct('implementation','C++','order',2));
T2 = fm(F,SPs,Dxyz,struct('implementation','C++','order',1));


%% Calculate exact answer
% Calculate grid node positions
[xx yy zz] = fmMeshGrid([n m o],Dxyz);

% Calculate one time/distance map for each SP
Texa = zeros(m,n,o,size(SPs,2));
for iter=1:size(SPs,2)
	Texa(:,:,:,iter) = sqrt((xx-SPs(1,iter)).^2 + ...
		(yy-SPs(2,iter)).^2 + (zz-SPs(3,iter)).^2);
end
Texa = min(Texa,[],4);


%% Calculate error
rel_errors = (T1-Texa)./(Texa);
per_pixel_rel_error = sum(abs(rel_errors(:)))/(m*n);

display(per_pixel_rel_error);

%% Plots
% Get non-singleton dimension
T1 = squeeze(T1);
T2 = squeeze(T2);
Texa = squeeze(Texa);
rel_errors = squeeze(rel_errors);

figure(1);

subplot(2,2,1);
imshow(T1,[],'Init','fit','XD',[1 Ly],'YD',[1 Lz]);
title('Fast marching solution');

subplot(2,2,2);
imshow(T2,[],'Init','fit','XD',[1 Ly],'YD',[1 Lz]);
title('Fast marching solution 2');

subplot(2,2,3);
imshow(Texa,[],'Init','fit','XD',[1 Ly],'YD',[1 Lz]);
title('Exact solution');

colorbar('Position', [.01 .1 .05 .8]);

subplot(2,2,4);
imshow(rel_errors,[],'Init','fit','XD',[1 Ly],'YD',[1 Lz]);
title('Relative errors');

colorbar('Position', [.9 .1 .05 .8]);

colormap(hot());
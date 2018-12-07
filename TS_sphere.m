% Sphere test for 3D fm. Shows cross sections

%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

% Domain dimensions
Lx = 8; 
Ly = 8;
Lz = 8;

% Number of grid nodes in each direction
n = Lx;
m = Ly;
o = Lz;

% dx, dy, dz
dx = Lx/n;
dy = Ly/m;
dz = Lz/o;
Dxyz = [dx dy dz];

% Speed map
F = ones(m,n,o);
% F = ones(m,n,o) + 4*rand(m,n,o);

% Source points
SPs = [Lx/2 Ly/2 Lz/2]';


%% Calculate T
T = fm(F,SPs, Dxyz, struct('implementation','Matlab','order',2));


%% Calculate exact answer
% Only works for unit F!!!
[xx yy zz] = fmMeshGrid([n m o], Dxyz);

% Calculate one time/distance map for each SP
Texa = zeros(m,n,o,size(SPs,2));
for iter=1:size(SPs,2)
	Texa(:,:,:,iter) = sqrt((xx-SPs(1,iter)).^2 + ...
		(yy-SPs(2,iter)).^2 + (zz-SPs(3,iter)).^2);
end
Texa = min(Texa,[],4);


%% Calculate error
rel_errors = (T-Texa)./(Texa);
per_pixel_rel_error = sum(abs(rel_errors(:)))/(m*n*o);

display(per_pixel_rel_error);


%% Plot
nCrossSections = 9;
CSinds = ceil(o./(nCrossSections:-1:1));

% calculate subplot numbers
xSubPlots = ceil(sqrt(nCrossSections));
ySubPlots = ceil(nCrossSections/xSubPlots);

% Get max values to enforce coherent colormaps
minT = min([T(:)' Texa(:)']);
maxT = max([T(:)' Texa(:)']);
minD = min(rel_errors(:));
maxD = max(rel_errors(:));

% Could also use these, instead of minD, maxD
mu = mean(rel_errors(:));
sig= std(rel_errors(:));

% Plot
for iter = 1:nCrossSections
	CS = CSinds(iter);
	% Plot cross section
	figure(1);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(T(:,:,CS),[minT maxT], ...
		'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
	figure(2);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(Texa(:,:,CS),[minT maxT], ...
		'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
	% Plot error cross section
	figure(3);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(rel_errors(:,:,CS),[minD maxD], ...
		'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
end

% Configure plots
figure(1); set(figure(1),'Name', 'FM CrossSections');
colormap(hot(300)); 
colorbar('Position', [.9 .1 .05 .8])
figure(2); set(figure(2),'Name', 'Exact Solution CrossSections');
colormap(hot(300));
colorbar('Position', [.9 .1 .05 .8])
figure(3); set(figure(3),'Name', 'Error CrossSections');
colormap(hot(300));
colorbar('Position', [.9 .1 .05 .8])
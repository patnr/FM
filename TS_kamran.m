%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

% Load Speed map
% load('data/kamran.mat','F')
load('data/Data.mat','F')

% Nondef nodes
F(F==-9999) = 0;

% Rearrange dimensions of F to get ordering (y, x, z).
% Use the cross section plots to determine which dimension should
% actually be the xyz dimensions, if this is not correct (ref: your
% email about C/Matlab dimension uncertainty).
% Note that whilst F is ordered (y,x,z), the vectors Nxyz, Lxyz, Dxyz
% and SPs should all be ordered (x,y,z).
F = permute(F,[2 3 1]);

% F = ones(size(F));
% F = ones([60 10 10]);

% Get dimensionality
Nyxz = size(F);
m = Nyxz(1); % y size
n = Nyxz(2); % x size
o = Nyxz(3); % z size
Nxyz = [n m o];
clear Nyxz;

% Lengths (not number of nodes) of dimensions
Lx = 100;
Ly = 100;
Lz = 100;

% Or: 
% Lx = n;
% Ly = m;
% Lz = o;

Lxyz = [Lx Ly Lz];

% Dxyz
Dxyz = Lxyz./Nxyz;

% Set source point
SPs = [Lx/2 Ly/2 Lz*(15/o)]';


%% Calculate T
T = fm(F,SPs, Dxyz, 'imp','C++','order',2);

% load('data/kamran.mat','T')


%% Plots
PlotF = 1;

% Calculate cross section indices for the different dimensions
nCrossSections = 9;
CSinds = ceil([m;n;o]*(1./(nCrossSections:-1:1)));

% calculate subplot numbers
xSubPlots = ceil(sqrt(nCrossSections));
ySubPlots = ceil(nCrossSections/xSubPlots);

minF = min(F(:));
maxF = max(F(:));

% Plot
for iter = 1:nCrossSections
	if PlotF
	% Plot cross sections of F -- Constant dim 1
	figure(1);
	CS = CSinds(1,iter);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(squeeze(F(CS,:,:)),[minF maxF], ...
		'Init','fit','XD',[1 Ly],'YD',[1 Lz]);
	
	% Plot cross sections of F -- Constant dim 2
	figure(2);
	CS = CSinds(2,iter);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(squeeze(F(:,CS,:)),[minF maxF], ...
		'Init','fit','XD',[1 Lx],'YD',[1 Lz]);
	
	% Plot cross sections of F -- Constant dim 3
	figure(3);
	CS = CSinds(3,iter);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(squeeze(F(:,:,CS)),[minF maxF], ...
		'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
	end
	
	% Plot cross sections of T
	figure(4);
	CS = CSinds(1,iter);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(squeeze(T(CS,:,:)),[], ...
		'Init','fit','XD',[1 Ly],'YD',[1 Lz]);
	
	figure(5);
	CS = CSinds(2,iter);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(squeeze(T(:,CS,:)),[], ...
		'Init','fit','XD',[1 Lx],'YD',[1 Lz]);
	
	figure(6);
	CS = CSinds(3,iter);
	subplot(ySubPlots, xSubPlots, iter)
	imshow(squeeze(T(:,:,CS)),[], ...
		'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
end

%% Configure plots
if PlotF
figure(1); set(figure(1),'Name', 'F CrossSections - y and z');
colormap(hot(300)); 
colorbar('Position', [.9 .1 .05 .8])
figure(2); set(figure(2),'Name', 'F CrossSections - x and z');
colormap(hot(300));
colorbar('Position', [.9 .1 .05 .8])
figure(3); set(figure(3),'Name', 'F CrossSections - x and y');
colormap(hot(300));
colorbar('Position', [.9 .1 .05 .8])
end


figure(4); set(figure(4),'Name', 'T CrossSections - y and z');
colormap(hot(300)); 
% colorbar('Position', [.9 .1 .05 .8]) for some reason colorbar
% doesn't work here.
figure(5); set(figure(5),'Name', 'T CrossSections - x and z');
colormap(hot(300)); 
% colorbar('Position', [.9 .1 .05 .8])
figure(6); set(figure(6),'Name', 'T CrossSections - x and y');
colormap(hot(300)); 
% colorbar('Position', [.9 .1 .05 .8])
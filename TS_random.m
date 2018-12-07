% Try out random speed maps, testing the stability of the
% method/implementation
% Also use a wall with speed map values = 0, in order to test
% stability.
% Try different amplitudes for the random noise, to get different
% warning messages.

%% Setup
% clear all;
close all;

% Add path of functions
addpath('./functions');


% Random with an obstacle
m  = 100; dy = 1; Ly = m*dy;
n  = 100; dx = 1; Lx = n*dx;
Dxyz = [dx dy];

% Random amplitude. Vary this relative to 1 in order to see the effect
% of extreme speed map values
RA = 3;

% Wall configuration
w = 2; l = floor(sqrt(m*n)/3); corner = floor([3/5*n, 3/5*m]);

F2 = ones(m,n) + RA*rand(m,n); 
F2(corner(2):corner(2)+l,corner(1):corner(1)+w) = 0;
F2(corner(2):corner(2)+w,corner(1):corner(1)+l) = 0;
% save('data/SpeedMaps.mat','F2')
% load('data/SpeedMaps.mat','F2');


%% Calc
SPs = [corner(1)+5, corner(2)+5]';
% SPs = getPoints(getNodes(SPs,Dxyz),Dxyz);

[T2 eFlag] = fm(F2,SPs,Dxyz,struct('implementation','C++','order',2));
% T2 = msfm2d(F2,SPs,1,0);


%% Plot
figure(1);
imshow(T2,[],'Init','fit','XD',[1 Lx],'YD',[1 Ly]);
title('Speed map with noise and a wall');
colormap(hot());
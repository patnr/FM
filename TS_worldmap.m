% Use a world map to create a speed map. 
% Set travel time on the ocean to something, and 0 on land. 
% Play around with opening/closing the different channels as well as
% different cities.

%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

% Load world map
World = im2double(imread('data/world.jpg'));
F3 = sum(World,3);

% Make ocean crossable and land uncrossable
xable = 1;
F3(F3<1.5) = 0; 
F3(F3>=1.5) = xable; 

% Extend world to *fully* (it needs to be more than one pixel wide if
% the rayTracing should succeed) open the Bering strait.
% F3 = [F3 xable*ones(size(F3,1),2)];

% Or close the Bering strait (uncomment the "open Bering strait" line)
F3 = F3(:,1:end-2);

% Open the Suez canal
xstart = 374;
ystart = 134;
for i=1:5
	F3(ystart,xstart:xstart+3) = xable;
	xstart = xstart + 1;
	ystart = ystart + 1;
end

% Open the panama canal
F3(172:172+2,172:172+2) = xable;

% Open the English channel
F3(97:97+1,320:320+1) = xable;

% Coordinates of
Oslo  = [335; 83];
Tokyo = [575; 120];
LA    = [104; 129];


%% Calculations
T3 = fm(F3,Oslo,[1 1],struct('implementation','C++','order',2));
% T3 = fm(F3,Oslo,[1 1],struct('implementation','mat','order',2));
RTsettings = struct('useFD',0,'stepSize',1);
path2 = rayTrace(T3,Tokyo,Oslo,RTsettings);


%% Plot
figure(1);

% Can use imshow
% imshow(T3,[],'Init','fit');

% Or contour. In this case, tweak T3 a little
T3(T3==Inf) = 1.05*max(T3(T3(:)~=Inf));
contourf(flipud(T3));

% Plot ray
hold on;
% If using contour, it's necessary to reflect y values
path2(1,:) = 305 + 1 - path2(1,:);
plot(path2(2,:),path2(1,:),'white','LineWidth',2);
hold off;

% Configure plot
cmp = colormap('jet');
cmp (end,:) = 0.3;
colormap(cmp);
title(['Distances from Oslo by boat' 10 'Shortest path to Tokyo']);
set(findall(gcf,'Type','text'),'FontSize',18);
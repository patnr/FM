% Solve a maze

%% Setup
clear all;
close all;

% Add path of functions
addpath('./functions');

F1 = 1000*im2double(imread('data/maze.gif'));
SPs = [10;10];
mazeEnd = [803; 803];


%% Calculations
T1 = fm(F1,SPs,[1 1],struct('implementation','C++','order',1));
settings = struct('useFD',0,'stepSize',2);
path1 = rayTrace(T1,mazeEnd,SPs,settings);


%% Plot
% Plot fast marching solution
figure(1);
imshow(T1,[0 1.1*max(T1(T1(:)~=Inf))],'Init','fit');

% Plot path
hold on;
plot(path1(2,:),path1(1,:),'white','LineWidth',2);
hold off;

% Colormap config
colormap('hot');
cmp=colormap();
cmp(end,:) = [0 0 1];
colormap(cmp);

% Title
title('Maze solver with ray trace');
set(findall(gcf,'Type','text'),'FontSize',18);
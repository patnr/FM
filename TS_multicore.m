% Test multicore speed gains

% Remember to start a new (or N-1 new, if you have N cores) Matlab
% processes (with -nojvm, if you like) with this directory as the PWD and
% addpath('../Multicore'), addpath('functions'), before running 
% startmulticoreslave('mcDir'). All this before running this script.
%
% Remember to set the nCores parameter in 'fm.m'.

%% Setup
clear all;
close all;

% Domain dimensions
Lx = 100; 
Ly = Lx;
Lz = Lx;

% Number of grid nodes in each direction
n = 20;
m = n;
o = n;

% Dxyz
dx = Lx/n;
dy = Ly/m;
dz = Lz/o;
Dxyz = [dx dy dz];

% Speed map
F = ones(m,n,o) + 0.2*rand(m,n,o);

% Source points
nSPs = 4;
SPs = diag([Lx Ly Lz]) * rand(3,nSPs);


%% Calculate Fast-marching solutions
% Run fm with one processor only
tic;
T1 = cell(1,nSPs);
eFlags1 = cell(1,nSPs);
for k=1:nSPs
	[T1{k} eFlags1{k}] = ...
		fm(F,SPs(:,k),Dxyz,'impl','C++','ord',1);
end
timeT1 = toc;

% Run fm on multiple processors
tic;
[T2 eFlags2] = fm(F,SPs,Dxyz,'impl','C++','useMC',1,'ord',1);
timeT2 = toc;

display([10 'Single time: ' num2str(timeT1)]);
display(['Multi  time: ' num2str(timeT2)]);
display(['Single/Multi time ratio: ' num2str(timeT1/timeT2) '.']);

display(['Discrepancy:' num2str(sum(T1{1}(:)-T2{1}(:)))]);
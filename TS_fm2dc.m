% Time comparison C++ code vs Matlab code

%% Setup
% Allow for non-semicolon-ended output
%#ok<*NOPTS>

% Add path of functions
addpath('./functions');

m = 40; % Num of y nodes
n = m;   % Num of x nodes

% Grid distances
dx = 1;
dy = 3;

% Speed map
F = ones(m,n);

% Source points
SPs = [4 4]';


%% Solve for T (distance map)
% Matlab version (with class heap implementation)
tic; 
% T1 = fm(F,SPs,[dx dy],'imp','mat','order',1); 
T1 = fm2d(F,SPs,dx,dy,int32(1)); 
T1time = toc 

% Matlab version (withOUT class heap implementation)
% NB: Outdated. Only implemented for order 1.
tic; 
T2 = fm2d_noHeap(F,SPs,dx,dy); 
T2time = toc

% C++ version
tic; 
% T3 = fm(F,SPs,[dx dy],'imp','C++','order',1);
T3 = fm2dc(F,SPs,dx,dy,int32(1)); 
T3time = toc

%% Ensure all answers are the same
T12discrepancy = sum(T1(:) - T2(:))
T13discrepancy = sum(T1(:) - T3(:))



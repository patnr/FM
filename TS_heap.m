% Test heap functionality

clear all;

% Add path of functions
addpath('./functions');

% Heap elements
N = 28;

% Create heap
myHeap = heap(N);
elements = N*rand(N,1);
for i=1:N
	myHeap.insert(elements(i),i);
end

% Print organized heap
myHeap.print;

% Swap and print disorganized
myHeap.swapElements(4,20);
myHeap.print;

% Reorganize, and print
myHeap.upHeap(20);
myHeap.print;
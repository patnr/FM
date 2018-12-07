% Compiles the c++ mex files
cd('functions');

% mex -c utils.cpp
% mex -c heap.cpp

mex fm2dc.cpp utils.cpp heap.cpp
mex fm3dc.cpp utils.cpp heap.cpp

cd('..');
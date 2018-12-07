function [T eFlag] = fm(F, SourcePoints, Dxyz, varargin)
% usage:
% [T eFlag] = fm(F, SourcePoints, Dxyz, opts)
%
%
% Implements the Fast-marching method to solve the Eikonal eqn
% in 2 or 3 dimensions:
% || grad(T) ||^2 = 1/F^2
%
% Note that this is an interface only. It treats the input and
% executes the correct functions accordingly.
%
%
% T				= distance (in time) map
% F				= velocity map
% SourcePoints	= matrix with dimensions (nDims * NumOfSPs)
% Lxyz			= Domain length in x y and z directions
%
% opts			= struct with the following case-sensitive fields:
%	'implementation' = 'Matlab', 'noHeap' or 'C++' 
%						'noHeap' is the old Matlab implementation, that 
%						doesn't use a heap to store narrow band pixels.
%	'useMC'		 = 1 or 0. Should the method use multicore?
%   'order'      = 1 or 2. Order accuracy of finite difference approx.
% opts can also be entered as name/value pairs
%
% The meaning of the error flag 'eFlag' is given by the displayed
% warnings (at the bottom of this file).


% Add path of functions
addpath('./functions');



if (nargin<3)
	error('Must supply at least the first 3 input arguments.');
end

% Set default options
opts.implementation = 'C++';
opts.useMC = 0;
opts.order = 2;
	
% Set input options
opts = parsepropval(opts,varargin{:});

% Check options for errors
opts.implementation = validatestring(opts.implementation,...
	{'Matlab','C++','noHeap'});
if ~any(opts.useMC == [0 1])
	error('''useMC'' must be either 0 or 1');
end
if ~any(opts.order == [1 2])
	error('''order'' must be either 1 or 2');
end

% Convert order to type integer, coz thats what the C++ program expects.
opts.order = int32(opts.order);



% Check dimensionality coherence
nDims = ndims(F);
if((nDims ~= 2) && (nDims ~= 3))
	error('Problem must be either 2- or 3-dimensional');
end
if((nDims ~= size(SourcePoints,1)) || ...
		(length(Dxyz) ~= nDims))
	error('Incompatible dimensionality of inputs');
end


% Get n, Lx, etc...
n = size(F,2);
m = size(F,1);
if (nDims==3)
	o = size(F,3);
end
dx = Dxyz(1);
dy = Dxyz(2);
if (nDims==3)
	dz = Dxyz(3);
end


% Check that F>=0
if any(F < 0)
	error('F must be >= 0');
end


% Check Nxyz and Lxyz conditions
if ((n < 1) || (m < 1) || (exist('o','var') && (o < 1)))
	error('n,m,o must all be > 0');
end
if ((dx <= 0) || (dy <= 0) || (exist('dz','var') && (dz <= 0)))
	error('Dxyz must all be > 0');
end


% Grid spacing
Lx = n*dx;
Ly = m*dy;
if (nDims==3)
	Lz = o*dz;
end


% Check that all input is real
if(~isreal(F) || ~isreal(SourcePoints) || ~isreal([dx dy]) || ...
		(exist('dz','var') && (~isreal(dz))))
	error('All inputs must be real');
end


% Check that SourcePoints are in the domain
if (	any(SourcePoints(:) < 0) || ...
		any(SourcePoints(1,:) > Lx) || ...
		any(SourcePoints(2,:) > Ly) || ...
		((nDims == 3) && any(SourcePoints(3,:) > Lz)))
	error('SourcePoints supplied are out of bounds');
end

% Make SourcePoints coherent with matlab dimension ordering.
% I.e. swap x and y rows (1 and 2). Convert to integer data type
SourcePoints(1:2,:) = flipud(SourcePoints(1:2,:));


% Run Multicore processes
if opts.useMC
	% Add multicore system path
	addpath('../Multicore');
	
	% Adjust nCores according to the number of cores in computer.
	nCores = 2;
	nSPs = size(SourcePoints,2);
	
	% make sure EvalTimeSingle is big enough...otherwise everything
	% will be done by one core. It has to be greater than the time
	% taken to calculate one time-distance map
	MCsettings.maxEvalTimeSingle = 200;
	
	% Should this be adjusted in case mod(nSPs,nCores) != 0 ???
	MCsettings.nrOfEvalsAtOnce = max(1,floor(nSPs/nCores));
	
	% Other settings
	MCsettings.multicoreDir = 'mcDir';
	MCsettings.useWaitbar = false;
	MCsettings.masterIsWorker = true;
	
	% Put correct parameters in the parameter cell
	parameterCell = cell(1,nSPs);
	if (nDims==2)
		switch opts.implementation
			case 'Matlab'
				funhandle = @fm2d;
			case 'C++'
				funhandle = @fm2dc;
			case 'noHeap'
				funhandle = @fm2d_noHeap;
		end
		for iSP=1:nSPs
			parameterCell{iSP} = ...
				{funhandle,F,SourcePoints(:,iSP),dy,dx,opts.order};
		end
	elseif (nDims==3)
		switch opts.implementation
			case 'Matlab'
				funhandle = @fm3d;
			case 'C++'
				funhandle = @fm3dc;
			case 'noHeap'
				error('noHeap has not been implemented in 3d');
		end
		for iSP=1:nSPs
			parameterCell{iSP} = ...
				{funhandle,F,SourcePoints(:,iSP),dy,dx,dz,opts.order};
		end
	end
	
	% Run multicore processes
	resultCells = startmulticoremaster(@encapsulateOutput, parameterCell, MCsettings);
	
	% Split result cells into T and eFlag cells. Create maxEFlag for
	% use with the warning messages given here.
	maxEFlag = 0;
	T = cell(nSPs,1);
	eFlag = cell(nSPs,1);
	for iSP=1:nSPs
		T{iSP} = resultCells{iSP}{1};
		eFlag{iSP} = resultCells{iSP}{2};
		if (eFlag{iSP} > maxEFlag)
			maxEFlag = eFlag{iSP};
		end
	end
	
% Run regular (non-multicore) processes
else
	if (nDims==2)
		switch opts.implementation
			case 'Matlab'
			[T eFlag] = fm2d(F,SourcePoints,dy,dx,opts.order);
			case 'C++'
			[T eFlag] = fm2dc(F,SourcePoints,dy,dx,opts.order);
			case 'noHeap'
				error('noHeap in 2D is outdated.');
% 			T = fm2d_noHeap(F,SourcePoints,dy,dx,opts.order);
		end
	elseif(nDims==3)
		switch opts.implementation
			case 'Matlab'
				[T eFlag] = fm3d(F,SourcePoints,dy,dx,dz,opts.order);
			case 'C++'
				[T eFlag] = fm3dc(F,SourcePoints,dy,dx,dz,opts.order);
			case 'noHeap'
				error('noHeap has not been implemented in 3D.');
		end
	end
	maxEFlag = eFlag;
end


% Display warning according to error flag.
if (maxEFlag==1)
	display([10 'Warning 1: Had to revert to order 1 scheme at least' 10 ...
		'once because the 2nd order scheme resultet in complex' 10 ...
		'values. Note that this does not mean the upwind' 10 ...
		'causality requirement was violated.']);
elseif (maxEFlag==2)
	display([10 'Warning 2: The upwind causality requirement has' 10 ...
		'been violated somehow. The algorithm had to use a' 10 ...
		'brute fix which might lead to large errors.' 10 ...
		'Too extreme speed map values input?']);
end

end
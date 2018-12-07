% Puts the output of the fast marching algorithms into a cell, so that
% it can be accepted and passed on by the multicore system (which only
% accepts one output argument).

function outCell = encapsulateOutput(fmXXX, F, SourcePoints, dy, dx, dz, order)
	if ((nargin<6) || (nargin>7))
		error('Wrong number of input arguments.');
	elseif (nargin==6)
		% 2D case
		order = dz;
		
		[T eFlag] = fmXXX(F,SourcePoints,dy,dx,order);
		outCell{1} = T;
		outCell{2} = eFlag;
		
	elseif (nargin==7)
		% 3D case
		[T eFlag] = fmXXX(F,SourcePoints,dy,dx,dz,order);
		outCell{1} = T;
		outCell{2} = eFlag;
	end
end
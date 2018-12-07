% Produces a grid whose nodes are centered on the pixels.
% Acts just as meshgrid.

function [xx yy zz] = fmMeshGrid(Nxyz,Dxyz)

nDims = length(Nxyz);
assert(nDims == length(Dxyz));
assert(nDims == nargout);

if (nDims==2)
	dx = Dxyz(1);
	dy = Dxyz(2);
	n  = Nxyz(1);
	m  = Nxyz(2);
	[xx yy] = meshgrid(dx*(1:n), dy*(1:m));
	xx = xx - dx/2;
	yy = yy - dy/2;
elseif (nDims==3)
	dx = Dxyz(1);
	dy = Dxyz(2);
	dz = Dxyz(3);
	n  = Nxyz(1);
	m  = Nxyz(2);
	o  = Nxyz(3);
	[xx yy zz] = meshgrid(dx*(1:n), dy*(1:m), dz*(1:o));
	xx = xx - dx/2;
	yy = yy - dy/2;
	zz = zz - dz/2;
else
	error('fmMeshGrid only supports 2D and 3D grids.');
end

end

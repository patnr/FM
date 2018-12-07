function Points = getPoints(Nodes, Dxyz)
% Does the opposite of getNodes, although getPoints applied to
% getNodes will generally not yield the original points, but the
% center of the pixel in which the original points were.
%
% Nodes  : nDims x nPoints
% Lengths :	1 x nDims
% nNodes  : 1 x nDims  (containing the number of nodes in each dim)

% Scale points
Points = diag(Dxyz) * Nodes;

% Subtract dx/2
for dim=1:length(Dxyz)
	Points(dim,:) = Points(dim,:) - Dxyz(dim)/2;
end

end
function Nodes = getNodes(Points, Dxyz)
% The nodes are situated in the center of its pixel.
% 
% The mapping, in one dimension, is as follows:
% if x is in ] (k-1)*dx ; k*dx ],    then x --> x_k
%
% The node positions  can easily be reconfigured by changing the
% initialization part of the fm algorithms. Also need to alter
% getPoints, which maps the opposite way (node --> coordinate) from 
% this function (pixel --> node). Whatever positioning choosen, it's
% important to be aware of it when comparing an fm solution to the
% exact solution.
%
% Points  : nDims x nPoints
% Lengths :	1 x nDims
% nNodes  : 1 x nDims  (containing the number of nodes in each dim)

Nodes = zeros(size(Points));
for dim = 1:size(Points,1);
	Nodes(dim,:) = ceil(Points(dim,:)/Dxyz(dim));
end

% Special case if any point coordinate == 0.
% Then set its node to 1.
Nodes(Points == 0) = 1;

end
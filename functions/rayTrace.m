function path = rayTrace(T, jiStart, jiEnd, settings)
% Traces a ray through a time-dist field from a starting point towards
% lower values.
%
% NB: Uses periodic BCs. => Requires walls around domain if you want
% the domain to be closed.
%
% Uses either the minimal neighbour or a finite difference scheme to
% calculate steps. If settings is used it should be a structure
% with two fields, 'useFD' = 1/0 and 'stepSize' = stepSize

%% Input treatment
if (nargin == 4)
	assert(isstruct(settings));
	assert(isfield(settings,'useFD'));
	assert(isfield(settings,'stepSize'));
	assert((settings.useFD == 1) || (settings.useFD == 0));
	assert(settings.stepSize > 0);
	
	useFD		= settings.useFD;
	stepSize	= settings.stepSize;
	
elseif (nargin == 3)
	if isstruct(jiEnd)
		error('Must supply and end (source) point');
	end
	useFD = 0;
	stepSize = 1;
else
	error('Must supply three or four input arguments.');
end

%% 
	m = size(T,1);
	n = size(T,2);
	
	% Establish initially allocated path length, and allocate
	% The factor 20 makes the init_path_length greater than the path
	% through the standard maze (=> avoids prompting user).
	numelT = numel(T);
	init_path_length = 20*ceil(sqrt(numelT));
	path_length = init_path_length;
	path = zeros(2,init_path_length);
	
	% Set starting point
	ijStart = jiStart(2:-1:1);
	path(:,1) = ijStart;
	newNode = getNodes(ijStart,[1 1]);
	
	% Get ending (source) point
	ijEnd = jiEnd(2:-1:1);
	
	
	% Offsets of neighbours
	offsets = [-1  1   0  0;
				0  0  -1  1];
	
	% Are we within stepSize of the end point?
	almostThere = 0;
	
	for iter=1:numelT-1
		% Reallocate or stop raytrace?
		if (iter>=path_length)
			display(['The path is now ' int2str(path_length) ...
				' pixels long. Do you wish to continue the ray tracing?']);
			yesno = input('(''y''/''n'') ');
			if ~any(strcmpi(yesno,{'y', 'yes'}))
				if ~almostThere
					display(['Warning: The ray tracing did not ' ...
						'converge to the end (source) point.']);
				end
				break;
			end				
			path = [path zeros(2,init_path_length)]; %#ok<AGROW>
			path_length = path_length + init_path_length;
		end
		
		% Get current pixel
		node = newNode;
		
		% Get neighbour coordinates
		neighs = repmat(node,1,4) + offsets;
		
		% Apply periodic BCs
		neighs = 1 + mod(neighs-1, repmat([m;n],1,4));
		
		% Convert to index
		neighs = sub2ind(neighs,m);
		
		% Take step
		if useFD
			% Calculate normalized gradient
			grad = [T(neighs(2))-T(neighs(1));  % Y/i downwards diff
					T(neighs(4))-T(neighs(3))]; % X/j rightwards diff
			
			% Although artificial, these two lines are necessary to
			% deal with infinite time-distances.
			grad(isinf(grad)) = sign(grad(isinf(grad))) * 1e100;
			grad(isnan(grad)) = 0; % Inf-Inf is set to zero gradient
			
			% Normalize
			grad = grad/norm(grad);
			
			% Calculate next pixel
			path(:,iter+1) = path(:,iter) - stepSize*grad;
			newNode = getNodes(path(:,iter+1),[1 1]);
		else
			% Choose minimal neighbour
			[~, minNeigh] = min(T(neighs(:)));
			
			% Add this neighbour to path
			% It is generally a bad idea to use stepSize ~= 1 when
			% using minimum method.
			path(:,iter+1) = path(:,iter) + stepSize*offsets(:,minNeigh);
			newNode = path(:,iter+1);
		end

		
		% Check for convergence
		if (almostThere || ((abs(ijEnd(1) - path(1,iter+1)) < stepSize) ...
				&& (abs(ijEnd(2) - path(2,iter+1)) < stepSize)))
			if almostThere
				if (T(sub2ind(newNode,m)) < T(sub2ind(node,m)))
					path = path(:,1:iter+1);
					return;
				else
					path = path(:,1:iter);
					return;
				end
			else
				almostThere = 1;
			end
		end
		
	end
end

% Customized conversions
function ind = sub2ind(ij, m)
	ind = zeros(1,size(ij,2));
	for p=1:size(ij,2)
		ind(p) = ij(1,p) + (ij(2,p)-1)*m;
	end
end
function ij = ind2sub(ind,m)
	ij = zeros(2,length(ind));
	for p=1:length(ind)
		ij(2,p) = ceil(ind/m);
		ij(1,p) = ind - (ij(2,p)-1)*m;
	end
end
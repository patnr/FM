function T = fm2d_noHeap(F,SourcePoints,dy,dx)

display(['NB: Outdated. Only implemented for order 1,' 10 ...
	'and doesn''t use the initialization' 10 ...
	'method of the other methods.']);
	
	
	%% Initialization
	
	% m is the y length, n is the x length (of domain)
	m = size(F,1);
	n = size(F,2);
	NumOfCPs = size(SourcePoints,2);
	
	% Time distance matrix. Initialize all points to -1
	T = zeros(m,n) - 1;
	
	% Matrix that keeps tabs on Frozen pixels
	Frozen = zeros(m,n);
	
	% Allocate heap.
	% The first row will contain the calculated travel times.
	% The 2nd row will contain the y and x indices (i,j)
	% reCPectively.
	Initial_Heap_Alloc_Size = ceil(m*n/10);
	Heap = zeros(2, Initial_Heap_Alloc_Size);
	
	% Initiate heap counter
	heapCount = 0;
	
	% Relative indices of the 4 neighbours to a given point
	iOffsets = [-1  1  0  0];
	jOffsets = [ 0  0 -1  1];
	
	% And then insert neighbours in narrow-band heap, calculating
	% distances. This is the initializing loop.
	for iter = 1:NumOfCPs
		% Freeze source points and calculate travel time
		CP = SourcePoints(:,iter);
		CPnode = getNodes(CP,[dy dx]);
		CPInd = sub2ind(CPnode(1),CPnode(2),m);
		T(CPInd) = norm(CP - getPoints(CPnode,[dy dx]))/F(CPInd);
		Frozen(CPInd) = 1;
		
		% For all neighbours of the source points
		for neigh = 1:4
			% Get index of neighbour. Store as i and j
			ni = CPnode(1) + iOffsets(neigh);
			nj = CPnode(2) + jOffsets(neigh);
			nInd = sub2ind(ni,nj,m);
			
			% If valid for consideration
			if (isInDomain(ni,nj,m,n) && ~Frozen(ni,nj))
				% If T(i,j) has not been previously calculated
				if (T(ni,nj) == -1)
					% Increase heap counter
					heapCount = heapCount + 1;
					if (heapCount > size(Heap,2))
						Heap = [Heap zeros(2,Initial_Heap_Alloc_Size)]; %#ok<AGROW>
					end
					% Add as last heap element
					heap_ind = heapCount;
					% Add pointer from T to the heap
					T(ni,nj) = heap_ind;
					% Add pointer from heap to T
					Heap(2,heap_ind) = sub2ind(ni,nj,m);
					% Calculate tentative T(neigh), 
					% and list it in the heap
					Heap(1,heap_ind) = norm(CP - getPoints([ni;nj],[dy dx]))/F(nInd);
				else
					% Do nothing. Because no pixels have been set to
					% "Frozen", no new information will be used in
					% calculating distances, and so the distance will
					% be the same as previously calculated.
				end
			end
		end
	end
	
	
	%% Loop
	% Now start the loop.
	% Loop until there are no more narrow band
	% neighbours, meaning the algorithm has finished.
	
	% RMQ: Could also use while any(~Frozen)  ???
	while (heapCount > 0)
		% Find min heap element and make it Frozen/Frozen. 
		% This involves registering it as Frozen and storing the time
		% in T instead of in the heap.
		[time heap_ind] = min(Heap(1,1:heapCount));
		CP = Heap(2,heap_ind);
		Frozen(CP) = 1;
		T(CP) = time;
		
		% Insert the last heap element where the CP was.
		% This keeps the heap array compact. Actually related
		% to heap-sort implementation stuff, although this method
		% is not used in Matlab.
		% NB: Don't do this if CP is the last heap element. 
		%     The pointer re-assign will then overwite the time stored
		%     in T(CP(1),CP(2)).
		if (heap_ind ~= heapCount)
			Heap(:,heap_ind) = Heap(:,heapCount);
			% Also remember the pointer in T to the heap element
			T(Heap(2,heap_ind)) = heap_ind;
		end
		
		% Decrement heap counter
		heapCount = heapCount - 1;

		
		% For all neighbours of CP
		for neigh = 1:4
			% Get neighbour coordinates
			[ni nj] = ind2sub(CP,m);
			ni = ni + iOffsets(1,neigh);
			nj = nj + jOffsets(1,neigh);
			
			% If (i,j) valid for consideration
			if (isInDomain(ni,nj,m,n) && ~Frozen(ni,nj))
				if (T(ni,nj) == -1)
					% Increase heap counter
					heapCount = heapCount + 1;
					if (heapCount > size(Heap,2))
						Heap = [Heap zeros(2,Initial_Heap_Alloc_Size)]; %#ok<AGROW>
					end
					
					% Add as last heap element
					heap_ind = heapCount;
					% Add pointer from T to the heap
					T(ni,nj) = heap_ind;
					% Add pointer from heap to T
					Heap(2,heap_ind) = sub2ind(ni,nj,m);
					% Calculate tentative T(neigh) distance,
					% and list it in the heap
					Heap(1,heap_ind) = ...
						calcDist(ni,nj,F(ni,nj),T,Frozen,m,n,dx,dy);
				else
					% Unlike the initialization loop, this time the
					% CP pixel is now new information that can
					% be used in distance calculation.
					
					% Get heap indice of (i,j)
					heap_ind = T(ni,nj);
					
					% Calculate distance
					Heap(1,heap_ind) = ...
						calcDist(ni,nj,F(ni,nj),T,Frozen,m,n,dx,dy);
				end
			end
		end
	end
end % fm2d


function dist = calcDist(i,j,Fij,T,Frozen,m,n,dx,dy)
	% As a tool, create a local 3x3 patch, where pixels are set to
	% Inf if they are outside the domain or not frozen.
	Patch = Inf(3);
	for pi=1:3
		for pj=1:3
			if (isInDomain(i-2+pi,j-2+pj,m,n) && Frozen(i-2+pi,j-2+pj))
				Patch(pi,pj) = T(i-2+pi,j-2+pj);
			end
		end
	end
	
	% Get the minimal vertical neighbour, or don't use vertical
	% component in gradient.
	if (isinf(Patch(1,2)) && isinf(Patch(3,2)))
		sy = 0;
		ymin = 0;
	else
		sy = 1;
		if (Patch(1,2) <= Patch(3,2))
			ymin = Patch(1,2);
		else
			ymin = Patch(3,2);
		end
	end
	% Get the minimal horizontal neighbour, or don't use vertical
	% component in gradient.
	if (isinf(Patch(2,1)) && isinf(Patch(2,3)))
		sx = 0;
		xmin = 0;
	else
		sx = 1;
		if (Patch(2,1) <= Patch(2,3))
			xmin = Patch(2,1);
		else
			xmin = Patch(2,3);
		end
	end
	
	% Calculate coefficients of 2nd degree polynomial correCPonding 
	% to the distance equation. Coeffs: aT^2 + bT + c = 0
	a = (sy/dy^2 + sx/dx^2);
	b = -2*(sy*ymin/dy^2 + sx*xmin/dx^2);
	c = (sy*ymin/dy)^2 + (sx*xmin/dx)^2 - (1/Fij)^2;
	
	% Solve quadratic equation. Only use the maximum root.
	dist = (-b+sqrt(b^2-4*a*c))/(2*a);
	
	% RMQ: by design of the algorithm, this root is always greater
	% than the neighbour pixel values, and always real. 
	% However, this I have only proven for the first-order difference
	% scheme. So if other schemes are used, you might want to do 
	% an error check here.
	% Error check example:
% 	if ~isreal(dist)
% 		diCPlay('Distance calculated from polynomial is complex');
% 	end	
end

function ind = sub2ind(i,j,m)
	ind = i + (j-1)*m;
end

function [i j] = ind2sub(ind,m)
	j = ceil(ind/m);
	i = ind - (j-1)*m;
end

function bool = isInDomain(i,j,m,n)
	bool = ( ( i >= 1) && (i <= m) ...
		&& (j >= 1) && (j <= n) );
end
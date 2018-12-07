function [T eFlag] = fm3d(F,SourcePoints, dy, dx, dz, order)
% See fm2d for better documentation, as the two scripts are pretty
% similar.

%% Initialization
m = size(F,1);
n = size(F,2);
o = size(F,3);
NumOfSPs = size(SourcePoints,2);

T = zeros(m,n,o) - 1;
Frozen = zeros(m,n,o);
narrowBand = heap(m*n*o);

% Relative indices of the 6 neighbours to a given point
iOffsets = [-1 1  0 0  0 0];
jOffsets = [ 0 0 -1 1  0 0];
kOffsets = [ 0 0  0 0 -1 1];

% Relative indices of the 26 neighbours to a given point
for i = -1:1
	for j = -1:1
		for k = -1:1
			ind = (k+2)+(j+1)*3+(i+1)*9;
			iOffsets_full(ind) = i; %#ok<AGROW>
			jOffsets_full(ind) = j; %#ok<AGROW>
			kOffsets_full(ind) = k; %#ok<AGROW>
		end
	end
end
iOffsets_full(14) = [];
jOffsets_full(14) = [];
kOffsets_full(14) = [];

% Error flag
eFlag = 0;

% Source point calculations
for iter = 1:NumOfSPs
	CP = SourcePoints(:,iter);
	
	CPnode = getNodes(CP,[dy dx dz]);
	CPInd = sub2ind(CPnode(1),CPnode(2),CPnode(3),m,n);
	
	T(CPInd) = norm(CP - getPoints(CPnode,[dy dx dz]))/F(CPInd);
	Frozen(CPInd) = 1;
	
	% In 3D we have 26 surrounding neighbours
	for neigh = 1:26
		ni = CPnode(1) + iOffsets_full(neigh);
		nj = CPnode(2) + jOffsets_full(neigh);
		nk = CPnode(3) + kOffsets_full(neigh);
		nInd = sub2ind(ni,nj,nk,m,n);
		
		if isInDomain(ni,nj,nk,m,n,o)
			time = norm(CP - getPoints([ni;nj;nk],[dy dx dz]))/F(nInd);
			if (T(nInd) >= 0)
				T(nInd) = min(time,T(nInd));
			else
				T(nInd) = time;
			end
			Frozen(nInd) = 1;
		end
	end
end

% Create initial narrow band
for ind = 1:m*n*o
	if Frozen(ind)
	
	[i j k] = ind2sub(ind,m,n);
		
	for neigh = 1:6
		ni = i + iOffsets(neigh);
		nj = j + jOffsets(neigh);
		nk = k + kOffsets(neigh);
		nInd = sub2ind(ni,nj,nk,m,n);
		
		if (isInDomain(ni,nj,nk,m,n,o) && ~Frozen(nInd))
			if (~narrowBand.isInHeap(nInd))
				[time tmpFlag] = calcTime(ni,nj,nk,F(nInd),T,Frozen,...
					m,n,o,dy,dx,dz,order);
				narrowBand.insert(time, nInd);
				eFlag = max(tmpFlag, eFlag);
			else
				% Do nothing.
			end
		end
	end
	end
end


%% Main Loop.
lCount = 0;
while (narrowBand.heapCount > 0)
	lCount = lCount + 1;
	
	[time CP] = narrowBand.getSmallest();
	[i j k] = ind2sub(CP,m,n);
	
	Frozen(CP) = 1;
	T(CP) = time;
	
	% For all neighbours of CP
	for neigh = 1:6
		% Get neighbour coordinates
		ni = i + iOffsets(neigh);
		nj = j + jOffsets(neigh);
		nk = k + kOffsets(neigh);
		nInd = sub2ind(ni,nj,nk,m,n);
		
		% If (i,j) valid for consideration
		if (isInDomain(ni,nj,nk,m,n,o) && ~Frozen(ni,nj,nk))
			if (~narrowBand.isInHeap(nInd))
				[time tmpFlag] = calcTime(ni,nj,nk,F(nInd),T,Frozen,...
					m,n,o,dy,dx,dz,order);
				narrowBand.insert(time, nInd);
			else
				[time tmpFlag] = calcTime(ni,nj,nk,F(nInd),T,Frozen,...
					m,n,o,dy,dx,dz,order);
				narrowBand.update(time, nInd);
			end
			eFlag = max(tmpFlag, eFlag);
		end
	end
end
end % fm2d



function [time eFlag] = calcTime(i,j,k,Fijk,T,Frozen,m,n,o,dy,dx,dz,order)

% Get patch
Patch = Inf([5 5 5]);
for pi=1:5
	for pj=1:5
		for pk=1:5
			if (isInDomain(i-3+pi,j-3+pj,k-3+pk,m,n,o) ...
					&& Frozen(i-3+pi,j-3+pj,k-3+pk))
				Patch(pi,pj,pk) = T(i-3+pi,j-3+pj,k-3+pk);
			end
		end
	end
end

% Inf short-circuit
if min([Patch(2,3,3) Patch(4,3,3) Patch(3,2,3) Patch(3,4,3) ...
		Patch(3,3,2) Patch(3,3,4)]) == Inf
	time = Inf;
	eFlag = 0;
	return;
end

% Min search
oy = 0;
ymin1 = 0;
ymin2 = 0;
if ~(isinf(Patch(2,3,3)) && isinf(Patch(4,3,3)))
	oy = 1;
	if (Patch(2,3,3) <= Patch(4,3,3))
		ymin1 = Patch(2,3,3);
		if ((order==2) && (Patch(1,3,3) <= Patch(2,3,3)))
			ymin2 = Patch(1,3,3);
			oy = 2;
		end
	else
		ymin1 = Patch(4,3,3);
		if ((order==2) && (Patch(5,3,3) <= Patch(4,3,3)))
			ymin2 = Patch(5,3,3);
			oy = 2;
		end
	end
end
% Min search
ox = 0;
xmin1 = 0;
xmin2 = 0;
if ~(isinf(Patch(3,2,3)) && isinf(Patch(3,4,3)))
	ox = 1;
	if (Patch(3,2,3) <= Patch(3,4,3))
		xmin1 = Patch(3,2,3);
		if ((order==2) && (Patch(3,1,3) <= Patch(3,2,3)))
			xmin2 = Patch(3,1,3);
			ox = 2;
		end
	else
		xmin1 = Patch(3,4,3);
		if ((order==2) && (Patch(3,5,3) <= Patch(3,4,3)))
			xmin2 = Patch(3,5,3);
			ox = 2;
		end
	end
end
% 3rd-dim min search
oz = 0;
zmin1 = 0;
zmin2 = 0;
if ~(isinf(Patch(3,3,2)) && isinf(Patch(3,3,4)))
	oz = 1;
	if (Patch(3,3,2) <= Patch(3,3,4))
		zmin1 = Patch(3,3,2);
		if ((order==2) && (Patch(3,3,1) <= Patch(3,3,2)))
			zmin2 = Patch(3,3,1);
			oz = 2;
		end
	else
		zmin1 = Patch(3,3,4);
		if ((order==2) && (Patch(3,3,5) <= Patch(3,3,4)))
			zmin2 = Patch(3,3,5);
			oz = 2;
		end
	end
end


% Coefficients of quadratic.
a = ((oy==1)+9/4*(oy==2))/dy^2 + ...
	((ox==1)+9/4*(ox==2))/dx^2 + ...
	((oz==1)+9/4*(oz==2))/dz^2;
b = ((-2*(oy==1)-6*(oy==2))*ymin1 + 3/2*(oy==2)*ymin2)/dy^2 + ...
	((-2*(ox==1)-6*(ox==2))*xmin1 + 3/2*(ox==2)*xmin2)/dx^2 + ...
	((-2*(oz==1)-6*(oz==2))*zmin1 + 3/2*(oz==2)*zmin2)/dz^2;
c = ((oy==1)*ymin1^2 + ...
	(oy==2)*(4*ymin1^2 + 1/4*ymin2^2 -2*ymin1*ymin2))/dy^2 + ...
	((ox==1)*xmin1^2 + ...
	(ox==2)*(4*xmin1^2 + 1/4*xmin2^2 -2*xmin1*xmin2))/dx^2 + ...
	((oz==1)*zmin1^2 + ...
	(oz==2)*(4*zmin1^2 + 1/4*zmin2^2 -2*zmin1*zmin2))/dz^2 - ...
	(1/Fijk)^2;

d = b^2-4*a*c;

% Error treatment and appropriate time-dist calculation
if ((d<0) && (order==2))
	order = 1;
	[time tmpFlag] = calcTime(i,j,k, Fijk, T, Frozen, m, n, o, dy, dx, dz, order);
	eFlag = max(1,tmpFlag);
elseif ((d<0) && (order==1))
	if (oy==0)
		ymin1 = Inf;
	end
	if (ox==0)
		xmin1 = Inf;
	end
	if (oz==0)
		zmin1 = Inf;
	end
	time = min([xmin1 ymin1 zmin1]) + 1/Fijk;
	eFlag = 2;
else
	% Solve quadratic equation. Only use the maximum root.
	time = (-b+sqrt(d))/(2*a);
	eFlag = 0;
end

end

% Index conversions
function ind = sub2ind(i,j,k,m,n)
ind = i + (j-1)*m + (k-1)*(m*n);
end

function [i j k] = ind2sub(ind,m,n)
k = ceil(ind/(m*n));
j = ceil((ind-(k-1)*(m*n))/m);
i = ind - (j-1)*m - (k-1)*(m*n);
end

function bool = isInDomain(i,j,k,m,n,o)
bool = ( ( i >= 1) && (i <= m) ...
	&& (j >= 1) && (j <= n) ...
	&& (k >= 1) && (k <= o) );
end
function [T eFlag] = fm2d(F,SourcePoints,dy,dx,order)

%% Initialization
% m is the y length, n is the x length (of domain)
m = size(F,1);
n = size(F,2);
NumOfSPs = size(SourcePoints,2);

% Time distance matrix. Initialize all points to -1
T = zeros(m,n) - 1;

% Matrix that keeps tabs on Frozen pixels
Frozen = zeros(m,n);

% Construct heap
narrowBand = heap(m*n);

% Relative indices of the 4 neighbours to a given point
iOffsets = [-1  1  0  0];
jOffsets = [ 0  0 -1  1];

% Relative indices of the 8 neighbours to a given point
iOffsets_full = [-1 -1 -1  0  0  1  1  1];
jOffsets_full = [-1  0  1 -1  1 -1  0  1];

% Initialize error flag
eFlag = 0;

% First, we calculate the time-distances to all 8 neighbouring pixels
% of all the source points. This calculation is done simply by taking
% the distance and dividing by the speed. We also freeze these.
% RMQ: Improve this by using the speed pixel of both CP and neighbour?
for iter = 1:NumOfSPs
	CP = SourcePoints(:,iter);
	
	% Get node and its indice
	CPnode = getNodes(CP,[dy dx]);
	CPInd = sub2ind(CPnode(1),CPnode(2),m);

	% Calculate travel time and freeze
	T(CPInd) = norm(CP - getPoints(CPnode,[dy dx]))/F(CPInd);
	Frozen(CPInd) = 1;
	
	% For all 8 neighbours of the source points
	for neigh = 1:8
		% Get index of neighbour. Store as i and j
		ni = CPnode(1) + iOffsets_full(neigh);
		nj = CPnode(2) + jOffsets_full(neigh);
		nInd = sub2ind(ni,nj,m);
		
		% Only check if isInDomain. (We'll recalculate frozen pixels)
		% If it happens to be another SP then that value should be lower
		% anyways. If it happens to be the already-calculated neighbour 
		% of another SP, then it's still fair game for a new calculation.
		if isInDomain(ni,nj,m,n)
			time = norm(CP - getPoints([ni;nj],[dy dx]))/F(nInd);
			if (T(nInd) >= 0)
				T(nInd) = min(time,T(nInd));
			else
				T(nInd) = time;
			end
			Frozen(nInd) = 1;
		end
	end
end

% Calculate the initial narrow band as all neighbouring pixels to the 
% ones that have been frozen. Note that this time, unlike the source-point
% loop, and henceforth in the algorithm, the neighbours of a pixel only 
% includes its 4 non-diagonal neighbours.
% Could do the loop without going through absolutely all pixels and
% the if-Frozen condition, but this would significantly convolute the code. 
% It's only O(N) costly.
for ind = 1:m*n
	if Frozen(ind)
	
	[i j] = ind2sub(ind,m);
		
	% For all 4 neighbours of the frozen points
	for neigh = 1:4
		% Get index of neighbour. Store as i and j
		ni = i + iOffsets(neigh);
		nj = j + jOffsets(neigh);
		nInd = sub2ind(ni,nj,m);
		
		% If (i,j) valid for consideration
		if (isInDomain(ni,nj,m,n) && ~Frozen(nInd))
			if (~narrowBand.isInHeap(nInd))
				[time tmpFlag] = calcTime(ni,nj,F(nInd),T,Frozen,m,n,dy,dx,order);
				narrowBand.insert(time, nInd);
				eFlag = max(tmpFlag, eFlag);
			else
				% Do nothing. As no new pixels have been frozen in
				% this loop, no new information is present, and hence
				% the value calculated will necessarily be the same.
			end
		end
	end
	end
end


%% Error check
% Check if complex values could arise as a consequence of the fact
% that the algorithm might be violated by the initialization (coming
% from extreme variations in the speed-map around the source points).
% Edit: No point in doing this, coz even with smooth F surrounding SPs
% the corner Frozen pixel may be furher away than a narrowBand
% pixel straight above, or to the right. Especially if the SP is
% not centered.


% if narrowBand.Times(1) < max(T(Frozen==1))
% % 	narrowBand.print()
% % 	max(T(Frozen==1)) 
% 	display(['Warning: After initialization, there are ' 10 ...
% 		'narrowBand pixels with smaller time-dist''s than ' 10 ...
% 		'some of the Frozen pixels. This is a violation of ' 10 ...
% 		'the upwind causality of the algorithm, and may ' 10 ...
% 		'yield invalid results. Consider smoothing the speed' 10 ...
% 		'map around the source point to avoid this.' 10]);
% end


%% Main Loop
% Now start the main loop.
% Loop until there are no more narrow band
% neighbours, meaning the algorithm has finished.
lCount = 0;
while (narrowBand.heapCount > 0)
	lCount = lCount + 1;

	% Get min heap element. 
	% This will be the new "center pixel" (CP).
	[time CP] = narrowBand.getSmallest();
	[i j] = ind2sub(CP,m);
	
	% Freeze and set time
	Frozen(CP) = 1;
	T(CP) = time;
	
	% For all neighbours of CP
	for neigh = 1:4
		% Get neighbour coordinates. 
		ni = i + iOffsets(1,neigh);
		nj = j + jOffsets(1,neigh);
		nInd = sub2ind(ni,nj,m);
		
		% If (i,j) valid for consideration
		if (isInDomain(ni,nj,m,n) && ~Frozen(nInd))
			if (~narrowBand.isInHeap(nInd))
				[time tmpFlag] = calcTime(ni,nj,F(ni,nj),T,Frozen,m,n,dy,dx,order);
				narrowBand.insert(time, nInd);
			else
				[time tmpFlag] = calcTime(ni,nj,F(ni,nj),T,Frozen,m,n,dy,dx,order);
				narrowBand.update(time,nInd);
			end
			
			eFlag = max(tmpFlag, eFlag);
			
%			Make a cool "video" showing the algorithm in action.
% 			showProgress(T,narrowBand, [i;j], [ni;nj], 0);
% 			pause(0.1);
% 			input('\n...');
% 			figure(1);
			
%			Debugging stuff
% 			unruly = narrowBand.checkHeapProperty(1);
% 			if unruly
% 				display('heapDisorganized');
% 				display(loopCount);
% 				input('\n...');
% 			end
% 			if isnan(time)
% 				display('Encountered nan.');
% 				lCount %#ok
% 				neigh  %#ok
% 				if input('\nEnter 0 to continue...'); return; end
% 			end
		end
	end
end
end % fm2d


%% calcTime
function [time eFlag] = calcTime(i,j,Fij,T,Frozen,m,n,dy,dx,order)
% Calculates the time-dist at (i,j) from the neighbouring Frozen
% pixels. The formula is the standard quadratic equation along with
% the min-max switches.


% To remove the cases [out of domain] and [non-frozen] from the
% rest of the code, create a local 3x3 patch, where pixels are set to
% Inf if it's any of these cases. Using Inf as a flag is useful,
% because it lets us use min() instead of tons of if-else's.
%
% Note that although it's unnecessary to establish the whole Patch,
% because only the cross (vertical and horizontal) pixels are used, it
% makes the code simpler. For speed, the C++ implementation should be
% used.
Patch = Inf(5);
for pi=1:5
	for pj=1:5
		if (isInDomain(i-3+pi,j-3+pj,m,n) && Frozen(i-3+pi,j-3+pj))
			Patch(pi,pj) = T(i-3+pi,j-3+pj);
		end
	end
end


% If all the surrounding cross values are Inf, then set time to Inf.
% Otherwise would result in d = b^2 - 4*a(=0) * c(=1/Fij=Inf) = NaN
if min([Patch(2,3) Patch(4,3) Patch(3,2) Patch(3,4)]) == Inf
	time = Inf;
	eFlag = 0;
	return;
end

% Get the minimal vertical neighbour, or don't use vertical
% component in gradient.
oy = 0;
ymin1 = 0;
ymin2 = 0;
if ~(isinf(Patch(2,3)) && isinf(Patch(4,3)))
	oy = 1;
	if (Patch(2,3) <= Patch(4,3))
		ymin1 = Patch(2,3);
		if ((order==2) && (Patch(1,3) <= Patch(2,3)))
			ymin2 = Patch(1,3);
			oy = 2;
		end
	else
		ymin1 = Patch(4,3);
		if ((order==2) && (Patch(5,3) <= Patch(4,3)))
			ymin2 = Patch(5,3);
			oy = 2;
		end
	end
end
% Get the minimal horizontal neighbour, or don't use vertical
% component in gradient.
ox = 0;
xmin1 = 0;
xmin2 = 0;
if ~(isinf(Patch(3,2)) && isinf(Patch(3,4)))
	ox = 1;
	if (Patch(3,2) <= Patch(3,4))
		xmin1 = Patch(3,2);
		if ((order==2) && (Patch(3,1) <= Patch(3,2)))
			xmin2 = Patch(3,1);
			ox = 2;
		end
	else
		xmin1 = Patch(3,4);
		if ((order==2) && (Patch(3,5) <= Patch(3,4)))
			xmin2 = Patch(3,5);
			ox = 2;
		end
	end
end


% Calculate coefficients of 2nd degree polynomial corresponding
% to the time-distance equation. Coeffs: aT^2 + bT + c = 0
% The (ox==1) and (ox==2) act as switches. Remove all the (ox==2) and
% what is left are the first order scheme coefficients.
a = ((oy==1)+9/4*(oy==2))/dy^2 + ((ox==1)+9/4*(ox==2))/dx^2;
b = ((-2*(oy==1)-6*(oy==2))*ymin1 + 3/2*(oy==2)*ymin2)/dy^2 + ...
	((-2*(ox==1)-6*(ox==2))*xmin1 + 3/2*(ox==2)*xmin2)/dx^2;
c = ((oy==1)*ymin1^2 + ...
	(oy==2)*(4*ymin1^2 + 1/4*ymin2^2 -2*ymin1*ymin2))/dy^2 + ...
	((ox==1)*xmin1^2 + ...
	(ox==2)*(4*xmin1^2 + 1/4*xmin2^2 -2*xmin1*xmin2))/dx^2 - ...
	(1/Fij)^2;

d = b^2-4*a*c;

% One can analytically prove that even if causality is not 
% violated (i.e. all Frozen pixels should always have lower values
% than all narrow band pixels), complex time-dist results may still 
% arise for 2nd order schemes. This is not the case for 1st order
% schemes, which may only yield complex results if the causality
% aspect of the algorithm has been violated. Higher order schemes
% are less tolerant of extreme variations in the speed-map values.
%
% This implementation first attempts to use the 1st order scheme if
% the second order scheme fails. If that also fails, meaning the 
% causality has been violated (which might happen if you're very 
% unlucky, because of previous reversions,or very extreme speed-map 
% differences around the source points) we simply add 1/Fij to the 
% smallest neighbour. In case of a reversion, an persistent error flag 
% 'eFlag' is set to 1. In case of a violation, it's set to 2. 
%
% This is better error treatment than what is used in the 'msfm'
% (multistencil fast-marching) code imho, where any complex
% calculation always results in [add 1/Fij to smallest neigh.], and
% which doesn't notify the user with any error flags.
if ((d<0) && (order==2))
	% Revert to order 1
	order = 1;
	[time tmpFlag] = calcTime(i,j, Fij, T, Frozen, m, n, dy, dx, order);
	eFlag = max(1,tmpFlag);
elseif ((d<0) && (order==1))
	% Add 1/Fij to smallest neighbour
	if (oy==0)
		ymin1 = Inf;
	end
	if (ox==0)
		xmin1 = Inf;
	end
	time = min(xmin1, ymin1) + 1/Fij;
	eFlag = 2;
else
	% All good.
	% Solve quadratic equation. Only use the maximum root.
	time = (-b+sqrt(d))/(2*a);
	eFlag = 0;
end

end

% These override the default functions, being a lot more efficient
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
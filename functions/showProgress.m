% Shows a pixel map containing both frozen and narrow band pixels for
% use with debugging. Search fm2d.m for "showProgress" to find an
% example for usage.

function Tshown = showProgress(T,narrowBand,CPij, Nij ,dontPlot)

	if (~exist('dontPlot','var') || nargin < 3)
		dontPlot = 0;
	end

	m = size(T,1);
	Tshown = T;
	nbx = zeros(2,narrowBand.heapCount);
	
	for k=1:narrowBand.heapCount
		ind = narrowBand.H2T(k);
		Tshown(ind) = narrowBand.Times(k);
		[i j] = ind2sub(ind,m);
		nbx(:,k) = [i;j];
	end
	
	if ~dontPlot
		imshow(Tshown,[-1 max(T(:))],'Init','fit');
		
		% Plot a cross in NB pixels
		for k=1:narrowBand.heapCount
			hold on;
			plot(nbx(2,:),nbx(1,:),'g*','MarkerSize',8);
			hold off;
		end
		
		% Plot neighbour pixel
		if exist('CPij','var')
			Nij = Nij(:);
			
			hold on;
			plot(Nij(2),Nij(1),'r*','MarkerSize',8);
			hold off;
		end
		
		% Plot current pixel
		if exist('CPij','var')
			CPij = CPij(:);
			
			hold on;
			plot(CPij(2),CPij(1),'b*','MarkerSize',8);
			hold off;
		end
	end

end


function [i j] = ind2sub(ind,m)
j = ceil(ind/m);
i = ind - (j-1)*m;
end

classdef heap < handle
	% See C++ code for explanations
	% NB: The formulae for indices for child/parent nodes
	% are slightly different in C++ since C++ arrays start at 0.
	
	properties
		% Arrays
		Times; % Holds tentative T values
		H2T;  % Pointers from the heap to T
		T2H;  % Pointers from T to the heap
		
		% Counters
		heapCount;
		Initial_Heap_Alloc_Size;
		allocatedSize;
	end
	
	methods
		% Constructor
		function obj = heap(N)
			obj.Initial_Heap_Alloc_Size = ceil(N/4);
			obj.heapCount = 0;
			obj.allocatedSize = obj.Initial_Heap_Alloc_Size;
			
			obj.Times = zeros(1,obj.Initial_Heap_Alloc_Size);
			obj.H2T  = zeros(1,obj.Initial_Heap_Alloc_Size);
			obj.T2H  = zeros(1,N);
		end		
	end
	
	methods
		
		function outInd = parentInd(~, inInd)
			outInd = floor(inInd/2);
		end
		
		function outInd = leftChildInd(~, inInd)
			outInd = 2*inInd;
		end
		
		function outInd = rightChildInd(~, inInd)
			outInd = 2*inInd + 1;
		end
		
		function outInd = lastParentInd(obj)
			outInd = floor(obj.heapCount/2);
		end
		
		
		
		function print(obj)
			% Constants
			S = 3; % Significant digits
			B = 4; % Whitespace between elements on last row
			
			nRows = 1+floor(log2(obj.heapCount));
			
			% Start recursive printing
			fprintf('\n');
			recurse(nRows,0,B);
			fprintf('\n\n');
			
			% Recursive print function
			function recurse(row,pad,spacing)
				if (row>1)
					newSpacing	= ceil(2*spacing + S);
% 					newSpacing	= ceil(2*spacing + S + 1);
					newPad		= ceil(pad + .5*spacing + .5*S);
					recurse(row-1,newPad,newSpacing);
				end
				
				padding    = repmat(' ',1,pad);
% 				whitespace = repmat(' ',1,spacing);
				
				fprintf('\n');
				fprintf(padding);
				for elem = 2^(row-1) : min((2^row-1),obj.heapCount)
					fprintf('%-*.*g', S+spacing, S, obj.Times(elem));
% 					fprintf(['%.' int2str(S) 'g' whitespace], ...
% 						obj.Times(elem));
				end
			end
		end
		
		function outInd = checkHeapProperty(obj, pInd)
		% Checks that the heap property is satisfied, from pInd and
		% down. Returns index in heap of the parent of the first
		% "unruly" child. Only returns the first such found.
		% Otherwise, returns 0.
		% Heap property: parent < child1, child2 for all parent nodes
			lChild = obj.leftChildInd(pInd);
			rChild = lChild + 1;
			
			outInd = 0;
			
			if((lChild <= obj.heapCount) && ...
					(obj.Times(lChild) < obj.Times(pInd)))
					outInd = pInd;
					return;
			end
			if((rChild <= obj.heapCount) && ...
					(obj.Times(rChild) < obj.Times(pInd)))
				outInd = pInd;
				return;
			end
			
			if ((outInd==0) &&  (lChild <= obj.lastParentInd()))
				outInd = obj.checkHeapProperty(lChild);
			end
			if ((outInd==0) && (rChild <= obj.lastParentInd()))
				outInd = obj.checkHeapProperty(rChild);
			end
		end

		
		
		function swapElements(obj,Ind1,Ind2)
			if(Ind1==Ind2)
				return;
			end
			
			% Swap Times
			tmp = obj.Times(Ind1);
			obj.Times(Ind1) = obj.Times(Ind2);
			obj.Times(Ind2) = tmp;
			
			% Swap T2H values
			% NB: Must come before H2T swaps
			obj.T2H(obj.H2T(Ind1)) = Ind2;
			obj.T2H(obj.H2T(Ind2)) = Ind1;
			
			% Swap H2T elems
			tmp = obj.H2T(Ind1);
			obj.H2T(Ind1) = obj.H2T(Ind2);
			obj.H2T(Ind2) = tmp;
		end
		
		function upHeap(obj, Ind)
			while(Ind>1)
				pInd = obj.parentInd(Ind);
				if (obj.Times(Ind) < obj.Times(pInd))
					obj.swapElements(Ind,pInd);
					Ind = pInd;
				else
					break;
				end
			end
		end
		
		function downHeap(obj, Ind)
			if (obj.heapCount < 2)
				return;
			end
			
			while (Ind <= obj.lastParentInd())
				child1 = obj.leftChildInd(Ind);
				child2 = obj.rightChildInd(Ind);
				minChild = Ind;
				
				if (obj.Times(child1) < obj.Times(Ind))
					minChild = child1;
				end
				if( (child2 <= obj.heapCount) && ...
						(obj.Times(child2) < obj.Times(minChild)))
					minChild = child2;
				end
				
				if (minChild ~= Ind)
					obj.swapElements(Ind,minChild);
					Ind = minChild;
				else
					break;
				end
			end
		end
		
		
		
		function bool = isInHeap(obj, Ind)
			bool = (obj.T2H(Ind) > 0);
		end
		
		function insert(obj, time, Ind)
			
			if(obj.heapCount == (obj.allocatedSize-1))
				obj.Times = [obj.Times zeros(1,obj.Initial_Heap_Alloc_Size)];
				obj.H2T  = [obj.H2T  zeros(1,obj.Initial_Heap_Alloc_Size)];
			end
			
			obj.heapCount = obj.heapCount + 1;
			
			obj.Times(obj.heapCount) = time;
			obj.H2T(obj.heapCount)   = Ind;
			obj.T2H(Ind)			 = obj.heapCount;
			
			obj.upHeap(obj.heapCount);
		end
			
		function update(obj, time, Ind)
			obj.Times(obj.T2H(Ind)) = time;
			obj.upHeap(obj.T2H(Ind));
		end
		
		function [time Ind] = getSmallest(obj)
			time = obj.Times(1);
			Ind  = obj.H2T(1);
			
			obj.swapElements(1,obj.heapCount);
			
			obj.heapCount = obj.heapCount - 1;
			
			obj.downHeap(1);
		end
	end
	
end
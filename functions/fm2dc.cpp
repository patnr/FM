#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "heap.hpp"
#include "utils.hpp"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

// Implements fm2d in C++ code. See fm2d.m for better documentation 
// about the actual workings of the algorithm.

// NB: When using matrix indices (i,j), they run from 1 to (m,n).
// When using 1-dim array index 'ind', it runs from 0 to (m*n-1).


// Calculate time-distances
double calcTime(int i, int j, double Fij, double *T, bool *Frozen, 
				int m, int n, double dy, double dx,
			 	int lCount, int order, int* eFlag);

// Implementation of fm2d -- 2-dim fast marching in C++
void fm2dc(double* T, double *F, double *SourcePoints, 
				double dy, double dx, int m, int n, int NumOfSPs, 
				int order, int* eFlag);



// Matlab gateway function
void mexFunction( int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[]) {
	
		// Size measures
		mwSize m,n,NumOfSPs;

		// Flag to signal if the computations resulted in complex values.
		int *eFlag;

		// Order
		int order;

		// Input values and pointers to input data
		double dx, dy;
		double *F, *T, *SourcePoints;
	
		// Get pointers to input matrices
		F = mxGetPr(prhs[0]);
		SourcePoints = mxGetPr(prhs[1]);
		
		// Get dimensionality
		m = mxGetM(prhs[0]);
		n = mxGetN(prhs[0]);
		NumOfSPs = mxGetN(prhs[1]);

		// Get dx and dy
		dy = mxGetScalar(prhs[2]);
		dx = mxGetScalar(prhs[3]);

		// Get order
		order = *((int *) mxGetData(prhs[4]));

		// Create output matrix and its data pointer
		plhs[0] = mxCreateDoubleMatrix(m,n, mxREAL);
		T = mxGetPr(plhs[0]);

		// Create error flag
		plhs[1] = myCreateIntScalar();
		eFlag = (int*) mxGetData(plhs[1]);

		// Run fast marching	
		fm2dc(T,F,SourcePoints,dy,dx,m,n,NumOfSPs,order,eFlag);

		return;
}



// Actual algorithm
void fm2dc(double* T, double *F, double *SourcePoints, 
				double dy, double dx, int m, int n, int NumOfSPs, 
				int order, int* eFlag) {
			
		bool *Frozen;
		double time;

		int lCount = 0;
		int tmpFlag=0;
		*eFlag = 0;
	
		// x and y positions of center points, their corresponding 
		// nodes and neighbour nodes.
		double	CPx, CPy, CPjx, CPiy, njx, niy;

		// Index variables of center points and neighbours
		int CPi, CPj, CPInd, ni, nj, nInd;


		// Index offsets of cross-neighbours
		int iOffsets[] = {-1, 1, 0, 0};
		int jOffsets[] = {0, 0, -1, 1};

		// Index offsets of Patch 
		int iOffsets_full[] = {-1, -1, -1,  0,  0,  1,  1,  1};
		int jOffsets_full[] = {-1,  0,  1, -1,  1, -1,  0,  1};

		// Allocate Frozen
		Frozen = (bool *) mxCalloc(m*n, sizeof(bool));

		// Initialize Heap
		heap narrowBand = heap(m*n);
		
		// Directly calculate the time-dist for the pixels surrounding SPs.
		for(int iter=0; iter<NumOfSPs; iter++) {
				CPy = SourcePoints[2*iter];
				CPx = SourcePoints[2*iter+1];
				CPi = getNode(CPy,dy); 
				CPj = getNode(CPx,dx);
				CPiy = getPoint(CPi,dy);
				CPjx = getPoint(CPj,dx);
				CPInd = ij2ind(CPi,CPj,m);
			
				// Calculate time-distance of source point node	
				T[CPInd] = sqrt(pow(CPiy-CPy,2) + pow(CPjx-CPx,2))/F[CPInd];
				Frozen[CPInd] = true;
				
				// For all neighbours
				for(int neigh=0; neigh<8; neigh++) {
						ni = CPi + iOffsets_full[neigh];
						nj = CPj + jOffsets_full[neigh];
						niy = getPoint(ni,dy);
						njx = getPoint(nj,dx);
						nInd = ij2ind(ni,nj,m);

						if (isInDomain(ni,nj,m,n)) {
								time = sqrt(pow(niy-CPy,2) + pow(njx-CPx,2))/F[nInd];
								if (Frozen[nInd]) {
										T[nInd] = MIN(time,T[nInd]);
								}
								else {
										T[nInd] = time;
								}
								Frozen[nInd] = true;
						}
				}
		}
		

		// Create the initial narrow band
		for(int ind=0; ind<(m*n); ind++) {
				if (Frozen[ind]) {
						CPi = ind2i(ind,m);
						CPj = ind2j(ind,m);

						for(int neigh=0; neigh<4; neigh++) {
								ni = CPi + iOffsets[neigh];
								nj = CPj + jOffsets[neigh];
								nInd = ij2ind(ni,nj,m);

								if(isInDomain(ni,nj,m,n) && !Frozen[nInd]) {
										if(!narrowBand.isInHeap(nInd)) {
												time = calcTime(ni,nj,F[nInd],T,Frozen,
																m,n,dy,dx,lCount, order, &tmpFlag);
												narrowBand.insert(time,nInd);
												if (tmpFlag>(*eFlag)) {
														*eFlag = tmpFlag;
												}
										}
								}
						}
				}
		}


		// Loop
		while(narrowBand.nElems() > 0) {
				lCount++;

				// Get minimal time-distance and index of this narrow-band element
				time = narrowBand.getSmallest(&CPInd);
				CPi = ind2i(CPInd,m);
				CPj = ind2j(CPInd,m);

				// Freeze and set time
				Frozen[CPInd] = true;
				T[CPInd] = time;

				// For all neighbours
				for(int neigh=0; neigh<4; neigh++) {
						ni 	= CPi + iOffsets[neigh];
						nj 	= CPj + jOffsets[neigh];
						nInd = ij2ind(ni,nj,m);

						// If valid for consideration
						if(isInDomain(ni,nj,m,n) && !Frozen[nInd]) {
								// If T(ni,nj) has not been previously calculated
								if(!narrowBand.isInHeap(nInd)) {
										time = calcTime(ni,nj,F[nInd],T,Frozen,
														m,n,dy,dx,lCount, order, &tmpFlag);
										narrowBand.insert(time,nInd);
								}
								else {
										time = calcTime(ni,nj,F[nInd],T,Frozen,
														m,n,dy,dx,lCount, order, &tmpFlag);
										narrowBand.update(time, nInd);
								}

								if (tmpFlag>(*eFlag)) {
										*eFlag = tmpFlag;
								}
						}
				}
		}

		// De-allocate memory
		mxFree(Frozen);
}


// Time-dist calculation
double calcTime(int i, int j, double Fij, double *T, bool *Frozen, 
				int m, int n, double dy, double dx,
			 	int lCount, int order, int* eFlag) {

		// Get infinity. 
		// NB: Is this good practice? As in the matlab code, using Inf as 
		// a flag simplifies the code, so that's why I use this. Maybe this
		// doesn't work with some compilers?
		double Inf = 1.0/0.0;

		// Temporary error flag
		int tmpFlag = 0;

		// NB: Don't change the ordering without also changing the if-else 
		// statements below that set xmin, ymin. 		
		int iOffsets[] = {-2, -1, 1, 2,    0,  0, 0, 0};
		int jOffsets[] = { 0,  0, 0, 0,   -2, -1, 1, 2};

		// Replaces Patch (of the Matlab implementation). In the C++ code,
		// we don't use a full patch, but only the values 
		// horizontal/vertical to the central pixel.
		double CrossVals[8];

		// time calculated
		double time;

		// Indices of neighbouring pixels to (i,j)
		int ni, nj, nInd;

		// Variables to calculate gradient
		double xmin1, xmin2, ymin1, ymin2;
		int ox, oy;

		// Variables for quadratic equation
		double a, b, c, d;

		// Get values of surrounding pixels
		for(int iter=0; iter<8; iter++) {
				ni = i + iOffsets[iter];
				nj = j + jOffsets[iter];
				nInd = ij2ind(ni,nj,m);
				if(isInDomain(ni,nj,m,n) && Frozen[nInd]) {
						CrossVals[iter] = T[nInd];
				}
				else {
						CrossVals[iter] = Inf;
				}
		}


		// Inf short-cicuit
		if ((CrossVals[1]==Inf) && (CrossVals[2]==Inf) && 
						(CrossVals[5]==Inf) && (CrossVals[6]==Inf)) {
				*eFlag = 0;
				time = Inf;
				return time;
		}


		// Get ymin
		oy = 0;
		ymin1 = 0;
		ymin2 = 0;
		if (!((CrossVals[1] == Inf) && (CrossVals[2] == Inf))) {
				oy = 1;
				if (CrossVals[1]<CrossVals[2]) {
						ymin1 = CrossVals[1];
						if ((order==2) && (CrossVals[0]<=CrossVals[1])) {
								ymin2 = CrossVals[0];
								oy = 2;
						}
				}
				else {
						ymin1 = CrossVals[2];
						if ((order==2) && (CrossVals[3]<=CrossVals[2])) {
								ymin2 = CrossVals[3];
								oy = 2;
						}
				}
		}
		// Get xmin
		ox = 0;
		xmin1 = 0;
		xmin2 = 0;
		if (!((CrossVals[5] == Inf) && (CrossVals[6] == Inf))) {
				ox = 1;
				if (CrossVals[5]<CrossVals[6]) {
						xmin1 = CrossVals[5];
						if ((order==2) && (CrossVals[4]<=CrossVals[5])) {
								xmin2 = CrossVals[4];
								ox = 2;
						}
				}
				else {
						xmin1 = CrossVals[6];
						if ((order==2) && (CrossVals[7]<=CrossVals[6])) {
								xmin2 = CrossVals[7];
								ox = 2;
						}
				}
		}


		// Calculate quadratic equation coefficients. 
		a = ((oy==1) + (oy==2)*9.0/4.0)/(dy*dy) +
				((ox==1) + (ox==2)*9.0/4.0)/(dx*dx); 
		b = (((oy==1)*(-2.0) + (oy==2)*(-6.0))*ymin1 + 
						(oy==2)*1.5*ymin2)/(dy*dy) + 
				(((ox==1)*(-2.0) + (ox==2)*(-6.0))*xmin1 + 
						(ox==2)*1.5*xmin2)/(dx*dx);
		c = ((oy==1)*ymin1*ymin1 + (oy==2)*(4.0*ymin1*ymin1 + 
								0.25*ymin2*ymin2 - 2.0*ymin1*ymin2))/(dy*dy) +
				((ox==1)*xmin1*xmin1 + (ox==2)*(4.0*xmin1*xmin1 + 
								0.25*xmin2*xmin2 - 2.0*xmin1*xmin2))/(dx*dx) -
			 	1.0/(Fij*Fij);

		d = (b*b) - (4.0*a*c);


		// Error treatment and appropriate time-dist calculation
		if((d<0.0) && (order==2)) {
				order = 1;
				time = calcTime(i,j,Fij,T,Frozen,m,n,dy,dx,lCount,order,&tmpFlag);
				*eFlag = MAX(1,tmpFlag);
		}
		else if ((d<0.0) && (order==1)) {
				if (oy==0)
						ymin1=Inf;
				if (ox==0)
						xmin1=Inf;

				*eFlag = 2;
				time = MIN(xmin1,ymin1) + 1.0/Fij;
		}
		else{
				*eFlag = 0;
				time = ((-b + sqrt(d))/(2.0*a));
		}
		
		return time;
}

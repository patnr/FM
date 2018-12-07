#include "mex.h"
#include "matrix.h"
#include <math.h>
#include "heap.hpp"
#include "utils.hpp"

#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))


// NB: When using matrix indices (i,j), they run from 1 to (m,n).
// When using 1-dim array index 'ind', it runs from 0 to (m*n-1).


// Calculate time-distances
double calcTime(int i, int j, int k, double Fijk, double *T, bool *Frozen, 
				int m, int n, int o, double dy, double dx, double dz,
				int lCount, int order, int* eFlag);

// Implementation of fm2d -- 2-dim fast marching in C++
void fm3dc(double* T, double *F, double *SourcePoints, 
				double dy, double dx, double dz, int m, int n, int o, int NumOfSPs,
				int order, int* eFlag);



// Matlab gateway function
void mexFunction( int nlhs, mxArray *plhs[], 
				int nrhs, const mxArray *prhs[]) {
	
		// Size measures
		mwSize m,n,o,NumOfSPs;

		// Flag to signal if the computations resulted in complex values.
		int *eFlag;

		// Order
		int order;

		// Pointer to the dimensions array
		const int * pd;

		// Input values and pointers to input data
		double dx, dy, dz;
		double *F, *T, *SourcePoints;;
	
		// Get pointers to input matrices
		F = mxGetPr(prhs[0]);
		SourcePoints = mxGetPr(prhs[1]);

		// Get dimensionality
		pd = mxGetDimensions(prhs[0]);
		m = pd[0];
		n = pd[1];
		o = pd[2];
		NumOfSPs = mxGetN(prhs[1]);

		// Get dx and dy
		dy = mxGetScalar(prhs[2]);
		dx = mxGetScalar(prhs[3]);
		dz = mxGetScalar(prhs[4]);

		// Get order
		order = *((int *) mxGetData(prhs[5]));

		// Create output matrix and its data pointer
		plhs[0] = mxCreateNumericArray(3, pd, mxDOUBLE_CLASS, mxREAL);
		T = mxGetPr(plhs[0]);
		
		// Create error flag
		plhs[1] = myCreateIntScalar();
		eFlag = (int*) mxGetData(plhs[1]);

		// Run fast marching	
		fm3dc(T,F,SourcePoints,dy,dx,dz,m,n,o,NumOfSPs,order,eFlag);

		return;
}



void fm3dc(double* T, double *F, double *SourcePoints, 
				double dy, double dx, double dz, int m, int n, int o, int NumOfSPs,
				int order, int* eFlag) {

		bool *Frozen;
		double time;

		// Debug stuff
		int lCount = 0;
		int tmpFlag=0;
		*eFlag = 0;
		
		// x and y Position of center points, their corresponding 
		// nodes and neighbour nodes.
		double CPy, CPx, CPz, CPiy, CPjx, CPkz, niy, njx, nkz;

		// Index variables of center points and neighbours
		int CPi, CPj, CPk, CPInd, ni, nj, nk, nInd;

		// Index offsets to get cross-neighbours
		int iOffsets[] = {-1, 1,  0, 0,  0, 0};
		int jOffsets[] = { 0, 0, -1, 1,  0, 0};
		int kOffsets[] = { 0, 0,  0, 0, -1, 1};

		int iOffsets_full[26];
		int jOffsets_full[26];
		int kOffsets_full[26];


		for(int i=-1;i<=1;i++) {
		for(int j=-1;j<=1;j++) {
		for(int k=-1;k<=1;k++) {
				int ind = (k+1) + (j+1)*3 + (i+1)*9;
				if(ind>13) {
						ind--;
				}

				iOffsets_full[ind] = i;
				jOffsets_full[ind] = j;
				kOffsets_full[ind] = k;
		}}}


		// Initialize Frozen
		Frozen = (bool *) mxCalloc(m*n*o, sizeof(bool));

		// Initialize Heap
		heap narrowBand = heap(m*n*o);

		// Do the initial narrow band calculations
		for(int iter=0; iter<NumOfSPs; iter++) {
				CPy = SourcePoints[3*iter];
				CPx = SourcePoints[3*iter+1];
				CPz = SourcePoints[3*iter+2];
				CPi = getNode(CPy,dy);
				CPj = getNode(CPx,dx);
				CPk = getNode(CPz,dz);
				CPiy = getPoint(CPi,dy);
				CPjx = getPoint(CPj,dx);
				CPkz = getPoint(CPk,dz);
				CPInd = ijk2ind(CPi,CPj,CPk,m,n);

				// Calculate time-distance of source point node	
				T[CPInd] = sqrt(pow(CPiy-CPy,2) + pow(CPjx-CPx,2)
							 + pow(CPkz-CPz,2))/F[CPInd];
				Frozen[CPInd] = true;
				
				// For all neighbours
				for(int neigh=0; neigh<26; neigh++) {
						ni = CPi + iOffsets_full[neigh];
						nj = CPj + jOffsets_full[neigh];
						nk = CPk + kOffsets_full[neigh];
						niy = getPoint(ni,dy);
						njx = getPoint(nj,dx);
						nkz = getPoint(nk,dz);
						nInd = ijk2ind(ni,nj,nk,m,n);
						

						if(isInDomain(ni,nj,nk,m,n,o)) {
								time = sqrt(pow(niy-CPy,2) + pow(njx-CPx,2) +
												pow(nkz-CPz,2))/F[nInd];
								if(Frozen[nInd]) {
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
		for(int ind=0; ind<(m*n*o); ind++) {
				if (Frozen[ind]) {
						CPi = ind2i(ind,m,n);
						CPj = ind2j(ind,m,n);
						CPk = ind2k(ind,m,n);

						for(int neigh=0; neigh<6; neigh++) {
								ni = CPi + iOffsets[neigh];
								nj = CPj + jOffsets[neigh];
								nk = CPk + kOffsets[neigh];
								nInd = ijk2ind(ni,nj,nk,m,n);

								if(isInDomain(ni,nj,nk,m,n,o) && !Frozen[nInd]) {
										if(!narrowBand.isInHeap(nInd)) {
												time = calcTime(ni,nj,nk, F[nInd], T, Frozen,
																m,n,o, dy,dx,dz, lCount, order, &tmpFlag);
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
				CPi = ind2i(CPInd,m,n);
				CPj = ind2j(CPInd,m,n);
				CPk = ind2k(CPInd,m,n);


				// Freeze and set time
				Frozen[CPInd] = true;
				T[CPInd] = time;

				// For all neighbours
				for(int neigh=0; neigh<6; neigh++) {
						ni 	= CPi + iOffsets[neigh];
						nj 	= CPj + jOffsets[neigh];
						nk 	= CPk + kOffsets[neigh];
						nInd = ijk2ind(ni,nj,nk,m,n);

						// If valid for consideration
						if(isInDomain(ni,nj,nk,m,n,o) && !Frozen[nInd]) {
								// If T(ni,nj) has not been previously calculated
								if(!narrowBand.isInHeap(nInd)) {
										time = calcTime(ni,nj,nk,F[nInd],T,Frozen,
														m,n,o,dy,dx,dz, lCount, order, &tmpFlag);
										narrowBand.insert(time,nInd);
								}
								else {
										time = calcTime(ni,nj,nk,F[nInd],T,Frozen,
														m,n,o,dy,dx,dz, lCount, order, &tmpFlag);
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





double calcTime(int i, int j, int k, double Fijk, double *T, bool *Frozen, 
				int m, int n, int o, double dy, double dx, double dz,
				int lCount, int order, int* eFlag) {

		// Get infinity
		double Inf = 1.0/0.0;

		// Temporary error flag
		int tmpFlag = 0;

		// NB: Don't change these without also changing the if-else 
		// statements below that set xmin, ymin.
		int iOffsets[] = {-2, -1, 1, 2,    0,  0, 0, 0,    0,  0, 0, 0};
		int jOffsets[] = { 0,  0, 0, 0,   -2, -1, 1, 2,    0,  0, 0, 0};
		int kOffsets[] = { 0,  0, 0, 0,    0,  0, 0, 0,   -2, -1, 1, 2};

		// Replaces Patch (of the Matlab implementation)
		double CrossVals[12];

		// time calculated
		double time;

		// Indices of neighbouring pixels to (i,j)
		int ni, nj, nk, nInd;

		// Variables to calculate gradient
		double xmin1, xmin2, ymin1, ymin2, zmin1, zmin2;
		int ox, oy, oz;

		// Variables for quadratic equation
		double a, b, c, d;

		for(int iter=0; iter<12; iter++) {
				ni = i + iOffsets[iter];
				nj = j + jOffsets[iter];
				nk = k + kOffsets[iter];
				nInd = ijk2ind(ni,nj,nk,m,n);
				if(isInDomain(ni,nj,nk,m,n,o) && Frozen[nInd]) {
						CrossVals[iter] = T[nInd];
				}
				else {
						CrossVals[iter] = Inf;
				}
		}


		// Infinity short-circuit
		if ((CrossVals[1]==Inf) && (CrossVals[2]==Inf) && 
						(CrossVals[5]==Inf) && (CrossVals[6]==Inf) &&
						(CrossVals[9]==Inf) && (CrossVals[10]==Inf)) {
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
		// Get zmin
		oz = 0;
		zmin1 = 0;
		zmin2 = 0;
		if (!((CrossVals[9] == Inf) && (CrossVals[10] == Inf))) {
				oz = 1;
				if (CrossVals[9]<CrossVals[10]) {
						zmin1 = CrossVals[9];
						if ((order==2) && (CrossVals[8]<=CrossVals[9])) {
								zmin2 = CrossVals[8];
								oz = 2;
						}
				}
				else {
						zmin1 = CrossVals[10];
						if ((order==2) && (CrossVals[11]<=CrossVals[10])) {
								zmin2 = CrossVals[11];
								oz = 2;
						}
				}
		}


		// Calculate quadratic equation coefficients. 
		a = ((oy==1) + (oy==2)*9.0/4.0)/(dy*dy) +
				((ox==1) + (ox==2)*9.0/4.0)/(dx*dx) +
				((oz==1) + (oz==2)*9.0/4.0)/(dz*dz); 
		b = (((oy==1)*(-2.0) + (oy==2)*(-6.0))*ymin1 + 
						(oy==2)*1.5*ymin2)/(dy*dy) + 
				(((ox==1)*(-2.0) + (ox==2)*(-6.0))*xmin1 + 
						(ox==2)*1.5*xmin2)/(dx*dx) +
				(((oz==1)*(-2.0) + (oz==2)*(-6.0))*zmin1 + 
						(oz==2)*1.5*zmin2)/(dz*dz);
		c = ((oy==1)*ymin1*ymin1 + (oy==2)*(4.0*ymin1*ymin1 + 
								0.25*ymin2*ymin2 - 2.0*ymin1*ymin2))/(dy*dy) +
				((ox==1)*xmin1*xmin1 + (ox==2)*(4.0*xmin1*xmin1 + 
								0.25*xmin2*xmin2 - 2.0*xmin1*xmin2))/(dx*dx) +
				((oz==1)*zmin1*zmin1 + (oz==2)*(4.0*zmin1*zmin1 + 
								0.25*zmin2*zmin2 - 2.0*zmin1*zmin2))/(dz*dz) -
			 	1.0/(Fijk*Fijk);

		d = (b*b) - (4.0*a*c);

		if((d<0.0) && (order==2)) {
				order = 1;
				time = calcTime(i,j,k,Fijk,T,Frozen,m,n,o,dy,dx,dz,
								lCount,order,&tmpFlag);
				*eFlag = MAX(1,tmpFlag);
		}
		else if ((d<0.0) && (order==1)) {
				if (oy==0)
						ymin1=Inf;
				if (ox==0)
						xmin1=Inf;
				if (oz==0)
						zmin1=Inf;

				*eFlag = 2;
				time = MIN((MIN(xmin1,zmin1)),(ymin1)) + 1.0/Fijk;
		}
		else{
				*eFlag = 0;
				time = ((-b + sqrt(d))/(2.0*a));
		}
		
		return time;
}

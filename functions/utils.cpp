#include <math.h>
#include "mex.h"

// Checks if Matlab-style indices (i,j) are within (1..m, 1..n)
// NB: Note that this 2-dimensional checking is
// impossible with single array indices...
bool isInDomain(int i, int j, int m, int n) {
	return ((i >= 1) && (j >= 1) && (i <= m) && (j <= n));	
}

bool isInDomain(int i, int j, int k, int m, int n, int o) {
	return ((i >= 1) && (j >= 1) && (k >= 1) && 
					(i <= m) && (j <= n) && (k <= o));	
}

// Print utilities for different types of matrices
void printMatrix(bool *Mat, mwSize m, mwSize n) {
	for(int i=0; i<m; i++) {
		mexPrintf("\n");
		for(int j=0; j<n; j++)
			mexPrintf("%d  ", *(Mat + j*m + i));
	}
	mexPrintf("\n");
}
void printMatrix(double *Mat, mwSize m, mwSize n) {
	for(int i=0; i<m; i++) {
		mexPrintf("\n");
		for(int j=0; j<n; j++)
			mexPrintf("%.3g  ", *(Mat + j*m + i));
	}
	mexPrintf("\n");
}
void printMatrix(int *Mat, mwSize m, mwSize n) {
	for(int i=0; i<m; i++) {
		mexPrintf("\n");
		for(int j=0; j<n; j++)
			mexPrintf("%d  ", *(Mat + j*m + i));
	}
	mexPrintf("\n");
}

// For use in debugging, prints the elements surrounding ind.
// w specifies the width/2 of the printout.
void printSurr(double *Mat, int w, int ind, int m) {
		mexPrintf("\n");
		int printInd;
		for(int i=1-w; i<w; i++) {
				mexPrintf("\n");
				for(int j=1-w; j<w; j++) {
						printInd = ind + m*j + i;
						if (printInd>=0)
								mexPrintf("%-6.3g ", Mat[printInd]);
						else
								mexPrintf("    ");
				}
		}
}
void printSurr(int *Mat, int w, int ind, int m) {
		mexPrintf("\n");
		int printInd;
		for(int i=1-w; i<w; i++) {
				mexPrintf("\n");
				for(int j=1-w; j<w; j++) {
						printInd = ind + m*j + i;
						if (printInd>=0)
								mexPrintf("%d ", Mat[printInd]);
						else
								mexPrintf("  ");
				}
		}
}
void printSurr(bool *Mat, int w, int ind, int m) {
		mexPrintf("\n");
		int printInd;
		for(int i=1-w; i<w; i++) {
				mexPrintf("\n");
				for(int j=1-w; j<w; j++) {
						printInd = ind + m*j + i;
						if (printInd>=0)
								mexPrintf("%d ", Mat[printInd]);
						else
								mexPrintf("  ");
				}
		}
}


// Convert Matlab indices to array index, and other way around
mwSize ij2ind(int i, int j, mwSize m) {
	return ((i-1) + (j-1)*m);
}
mwSize ind2i(int ind, mwSize m) {
	return 1 + ind - (ind/m)*m;
}
mwSize ind2j(int ind, mwSize m) {
	return ((ind)/m) + 1;
}

mwSize ijk2ind(int i, int j, int k, mwSize m, mwSize n) {
		return ((k-1)*(m*n) + (j-1)*m + (i-1));
}
mwSize ind2i(int ind, mwSize m, mwSize n) {
		return 1 + ind - ((ind - (ind/(m*n))*(m*n))/m)*m - (ind/(m*n))*(m*n);
}
mwSize ind2j(int ind, mwSize m, mwSize n) {
		return 1 + (ind - (ind/(m*n))*(m*n))/m;
}
mwSize ind2k(int ind, mwSize m, mwSize n) {
		return 1 + ind/(m*n);
}



// Mappings between nodes and pixels and points
int getNode(double x, double dx) {
		int i = (int) ceil(x/dx);
		if (i==0)
				i = 1;
		return i;
}

double getPoint(int node, double dx) {
		return ((double) node)*dx - dx/2.0;
}

// Create mxArray for holding an integer scalar
mxArray * myCreateIntScalar() {
		int dimsOfAScalar[2] = {1,1};
		return mxCreateNumericArray(1,dimsOfAScalar,mxINT32_CLASS,mxREAL);
}

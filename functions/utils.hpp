#include "mex.h"

// Checks if Matlab-style indices (i,j) are within (1..m, 1..n)
// NB: Note that this 2-dimensional checking 2-dimensional checking is
// impossible with single array indices...
bool isInDomain(int i, int j, int m, int n);
bool isInDomain(int i, int j, int k, int m, int n, int o);

// Print utilities for different types of matrices
void printMatrix(double *Mat, mwSize m, mwSize n);
void printMatrix(int *Mat, mwSize m, mwSize n);
void printMatrix(bool *Mat, mwSize m, mwSize n);

// Print surrounding area
void printSurr(double *Mat, int w, int ind, int m);
void printSurr(int *Mat, int w, int ind, int m);
void printSurr(bool *Mat, int w, int ind, int m);

// Convert Matlab indices to array index, and other way around
mwSize ij2ind(int i, int j, mwSize m);
mwSize ind2i(int ind, mwSize m);
mwSize ind2j(int ind, mwSize m);

// 3D versions of the same -- distinguished by signature
mwSize ijk2ind(int i, int j, int k, mwSize m, mwSize n);
mwSize ind2i(int ind, mwSize m, mwSize n);
mwSize ind2j(int ind, mwSize m, mwSize n);
mwSize ind2k(int ind, mwSize m, mwSize n);

// Mappings between nodes and pixels and points
int getNode(double x, double dx);
double getPoint(int node, double dx);

// Create mxNumericArray scalar
mxArray * myCreateIntScalar();


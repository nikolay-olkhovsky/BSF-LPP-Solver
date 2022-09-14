#pragma once
#include "Problem-Parameters.h"		// Problem Parameters 

//=========================== User-defined Problem Types =========================
typedef double PT_float_T;
typedef PT_float_T PT_vector_T[PP_N];			// Vector of dimension n
typedef PT_float_T	PT_matrix_T[PP_MM][PP_N + 1];
typedef PT_float_T	PT_column_T[PP_MM];

//=========================== User-defined BSF Types =========================
struct PT_bsf_parameter_T {		// Type of Parameter for workers (current approximation)
	PT_vector_T x;		// Current approximation
	int indexToBlock;	// Index of variable to block
	int m;				// Index for new inequality (must be > 0)
	int sign;			// Sign of new inequality
	PT_float_T b;		// Right part of new inequality
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	PT_float_T* a;				// Pointer to row: a[0],...,a[n-1]
	PT_float_T* b;				// Pointer to constant term b: a_0x_0+...+a_{n-1}x_{n-1} \leq b
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T projection;	// Point projection onto hyperplane
	int nonZeroCounter;
	bool pointIn;
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	bool pointIn;
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	bool pointIn;
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
// Not used
};

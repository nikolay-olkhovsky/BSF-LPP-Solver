#pragma once
#include "CBSFBase.h"
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string>

//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)

using namespace std;

class CBSFTask : public CBSFBase
{
public:
	//========================== Problem variables ====================================
	int		PD_m;						// Current number of inequalities
	int		PD_n;						// Curtrent space dimension
	int		PD_numShiftsSameLength;		// Number of shifts with the same length
	int		PD_numDetDir;				// Number of sequential states "Determine Direction"
	bool		PD_pointIn;					// Point is inside polytope
	double	PD_shiftLength;				// Shift length
	int		PD_state;					// State of Job Dispatcher (see PC_bsf_JobDispatcher)
	int		PD_indexToBlock;			// Index of variable to block
	int		PD_firstZcvI;				// Index of first zero-costs variable
	int		PD_firstLcvI;				// Index of first low-costs variable
	double	PD_minCostPercentile;		// Minimum cost percentile
	double	PD_objF_u;
	double	PD_objF_w;
	////========================== Problem structures ====================================
	static PT_matrix_T PD_A;					// Matrix of coefficients of inequalities 
	static PT_column_T PD_b;					// Column of the constant terms of the system Ax <= PD_b
	static PT_vector_T PD_c;					// Objective Function Coefficients
	static PT_vector_T PD_apexPoint;			// Apex point
	static PT_vector_T PD_u;					// Base point on Polytope
	static PT_vector_T PD_direction;			// Unit vector to set shift direction
	static PT_vector_T PD_hi;					// Higher bound
	static PT_vector_T PD_lo;					// Lower bound
	static PT_vector_T PD_unitObjVector;		// = c/||c||
	static PT_vector_T PD_objVector;			// = PD_unitObjectiveVector * PP_OBJECTIVE_VECTOR_LENGTH
	int PD_objI[PP_N];					// Index of objective variables in absolute descending order
	//========================== Input/Output ====================================
	string PD_problemName;
	//
	//
	// 
	//========================== Matrix format ====================================
	// nor - number of rows; noc - number of columns; non - number of non-zero values
	string PD_MTX_File_A; /* File of matrix A (only non-zero elements):
	------------ begin of file -------------
	nor noc non // nor=m (number of inequalities); noc=n (dimension); non - number of non-zero values
	i_1 j_1 A_{{i_1}{j_1}} // i_1 - index of row; j_1 - index of column
	...
	i_k j_k A_{{i_k}{j_k}}
	------------ end of file----------------*/
	string PD_MTX_File_b; /* File of column b:
	------------ begin of file -------------
	nor noc // nor=m (number of inequalities); noc=1
	b_1
	...
	b_{nor}
	------------ end of file----------------*/
	string PD_MTX_File_lo; /* File of variable lower bounds:
	------------ begin of file -------------
	nor noc // nor=n (dimension); noc=1
	lo_1
	...
	lo_{nor}
	------------ end of file----------------*/
	string PD_MTX_File_hi; /* File of variable higher bounds:
	------------ begin of file -------------
	nor noc // nor=n (dimension); noc=1
	lo_1
	...
	lo_{nor}
	------------ end of file----------------*/
	string PD_MTX_File_c; /* File of variable higher bounds:
	------------ begin of file -------------
	nor noc // nor=n (dimension); noc=1
	c_1
	...
	c_{nor}
	------------ end of file----------------*/
	string PD_MTX_File_x0;/* File of starting point in the following format:
	------------ begin of file -------------
	nor noc // nor=n (dimension); noc=1
	x_1
	...
	x_{nor}
	------------ end of file----------------*/
	string PD_MTX_File_so;/* Solution file in the following format:
	------------ begin of file -------------
	nor noc // nor=n (dimension); noc=1
	x_1
	...
	x_{nor}
	------------ end of file----------------*/

	void PC_bsf_CopyParameter(
		PT_bsf_parameter_T parameterIn,
		PT_bsf_parameter_T* parameterOutP
	);
	void PC_bsf_Init(
		bool* success
	);
	void PC_bsf_IterOutput(
		PT_bsf_reduceElem_T* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	);
	void PC_bsf_IterOutput_1(
		PT_bsf_reduceElem_T_1* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	);
	void PC_bsf_IterOutput_2(
		PT_bsf_reduceElem_T_2* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	);
	void PC_bsf_IterOutput_3(
		PT_bsf_reduceElem_T_3* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	);
	void PC_bsf_JobDispatcher(
		PT_bsf_parameter_T* parameter,
		int* job,
		bool* exit,
		double t
	);
	//
	void PC_bsf_MapF(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T* reduceElem,
		int* success
	);
	void PC_bsf_MapF_1(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T_1* reduceElem,
		int* success
	);
	void PC_bsf_MapF_2(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T_2* reduceElem,
		int* success
	);
	void PC_bsf_MapF_3(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T_3* reduceElem,
		int* success
	);
	void PC_bsf_ParametersOutput(
		PT_bsf_parameter_T parameter
	);
	void PC_bsf_ProblemOutput(
		PT_bsf_reduceElem_T* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	);
	void PC_bsf_ProblemOutput_1(
		PT_bsf_reduceElem_T_1* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	);
	void PC_bsf_ProblemOutput_2(
		PT_bsf_reduceElem_T_2* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	);
	void PC_bsf_ProblemOutput_3(
		PT_bsf_reduceElem_T_3* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	);
	void PC_bsf_ProcessResults(
		PT_bsf_reduceElem_T* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	);
	void PC_bsf_ProcessResults_1(
		PT_bsf_reduceElem_T_1* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	);
	void PC_bsf_ProcessResults_2(
		PT_bsf_reduceElem_T_2* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	);
	void PC_bsf_ProcessResults_3(
		PT_bsf_reduceElem_T_3* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	);
	void PC_bsf_ReduceF(
		PT_bsf_reduceElem_T* x,
		PT_bsf_reduceElem_T* y,
		PT_bsf_reduceElem_T* z
	);
	void PC_bsf_ReduceF_1(
		PT_bsf_reduceElem_T_1* x,
		PT_bsf_reduceElem_T_1* y,
		PT_bsf_reduceElem_T_1* z
	);
	void PC_bsf_ReduceF_2(
		PT_bsf_reduceElem_T_2* x,
		PT_bsf_reduceElem_T_2* y,
		PT_bsf_reduceElem_T_2* z
	);
	void PC_bsf_ReduceF_3(
		PT_bsf_reduceElem_T_3* x,
		PT_bsf_reduceElem_T_3* y,
		PT_bsf_reduceElem_T_3* z
	);
	void PC_bsf_SetListSize(
		int* listSize
	);
	void PC_bsf_SetInitParameter(
		PT_bsf_parameter_T* parameter
	);
	void PC_bsf_SetMapListElem(
		PT_bsf_mapElem_T* elem,
		int i
	);
	//====================== Problem Functions ===========================
	bool		Conversion();
	void		DetermineDirection(PT_bsf_parameter_T* parameter, bool* exit, bool* repeat);
	PT_float_T	Distance(PT_vector_T x, PT_vector_T y);
	bool		GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector);
	bool		LoadMatrixFormat();
	void		MakeObjVector(PT_vector_T c, PT_vector_T objVector);
	void		MovingOnSurface(PT_vector_T ptr_unitVectorToSurface, PT_vector_T basePoint, PT_vector_T x, bool* goOn);
	PT_float_T	ObjF(PT_vector_T x);
	void		ObjUnitVector(PT_vector_T objUnitVector);
	//
	bool		PointInHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b);
	bool		PointInHalfspace_s(PT_vector_T point, PT_vector_T a, PT_float_T b);
	bool		PointInPolytope_s(PT_vector_T x);
	void		ProblemOutput(double elapsedTime);
	double		ProblemScale();
	bool		SavePoint(PT_vector_T x, const char* filename, double elapsedTime);
	void		Shift(PT_vector_T basePoint, PT_vector_T direction, PT_float_T siftLength, PT_vector_T endPoint);
	void		ShrinkUnitVector(PT_vector_T objUnitVector, int shrinkBound);
	void		SkipComments(FILE* stream);
	void		SortObjVarI();
	void		Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z);
	void		Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint);
	void		Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y);
	void		Vector_DivideEquals(PT_vector_T x, double r);
	PT_float_T	Vector_DotProductSquare(PT_vector_T x, PT_vector_T y);
	void		Vector_EpsZero(PT_vector_T x);
	PT_float_T	Vector_Norm(PT_vector_T x);
	PT_float_T	Vector_NormSquare(PT_vector_T x);
	void		Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector);
	void		Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y);
	void		Vector_MultiplyEquals(PT_vector_T x, double r);
	void		Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector);
	void		Vector_ProjectOnHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b, PT_vector_T projection, int* exitCode);
	void		Vector_ResetToZero(PT_vector_T x);
	void		Vector_Round(PT_vector_T x);
	void		Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
	void		Vector_Unit(PT_vector_T vector);
	//
};


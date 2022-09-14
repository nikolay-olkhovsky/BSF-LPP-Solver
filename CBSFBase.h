#pragma once
#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include "Problem-Parameters.h"		// Problem Parameters 
#include "Problem-Types.h"

// ========================= Constants ========================
#define BD_EXIT	true			// End of Work
//======================= Compatibility ==========================
#ifndef PP_BSF_MAX_MPI_SIZE
#ifdef PP_MAX_MPI_SIZE
#define PP_BSF_MAX_MPI_SIZE PP_MAX_MPI_SIZE
#endif // PP_MAX_MPI_SIZE
#endif // PP_BSF_MAX_MPI_SIZE


//=========================== BSF Types ===============================
struct BT_order_T {
	char exit;		// true, if worker must stop
	int jobCase;
	int iterCounter;
	PT_bsf_parameter_T parameter;
};

struct BT_extendedReduceElem_T {// Extended element type of reduce list
	PT_bsf_reduceElem_T elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};

struct BT_extendedReduceElem_T_1 {// Extended element type of reduce list
	PT_bsf_reduceElem_T_1 elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};

struct BT_extendedReduceElem_T_2 {// Extended element type of reduce list
	PT_bsf_reduceElem_T_2 elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};

struct BT_extendedReduceElem_T_3 {// Extended element type of reduce list
	PT_bsf_reduceElem_T_3 elem;	// Element of reduce list
	int reduceCounter;			// Reduce Counter
};

using namespace std;

class CBSFBase
{
public:
	/*---------------- BSF Skeleton Variables ----------------------*/
	int BSF_sv_addressOffset;			// Contains the number of the first element of the map-sublist scheduled to the current worker process.
	int BSF_sv_iterCounter;				// Contains the number of iterations performed so far.
	int BSF_sv_jobCase;					// Contains the number of the current activity (job) in workflow.
	int BSF_sv_mpiRank;					// Contains the rank of current MPI process.
	int BSF_sv_mpiMaster;				// Contains the rank of master MPI process.
	int BSF_sv_numberInSublist;			// Contains the relative number of the element in the map-sublist that the function Map is currently applied to. 
	int BSF_sv_numOfWorkers;				// Contains the total number of the worker processes.
	PT_bsf_parameter_T BSF_sv_parameter; // Contains the order parameters.
	int BSF_sv_sublistLength;			// Contains the length of a map-sublist scheduled to a worker process.
	// ========================= Variables ========================
	static int	BD_elemsPerWorker;	// Number of elements per worker
	static bool BD_exit;			// Indicator of ending the calculations
	static int	BD_jobCase;		// Default value
	static int	BD_listSize;		// Size of the list
	static int	BD_masterRank;		// Rank of Master = Number of MPI processes - 1
	static int	BD_newJobCase;	// Default value
	static int	BD_numOfWorkers;	// Number of Workers = Number of MPI processes - 1
	static int	BD_rank;			// Rank of the current process MPI 
	static int	BD_size;			// Number of MPI processes
	static bool BD_success;			// Indicator of Successfulness of Initialization
	static int	BD_tailLength;		// Length of the remainder of the list after dividing by the number of workers
	// ========================= Time Variables ========================
	static double BD_t;				// Total time
	static int BD_iterCounter;
	// ========================= Lists ========================
	static BT_extendedReduceElem_T* BD_extendedReduceList;
	static BT_extendedReduceElem_T_1* BD_extendedReduceList_1;
	static BT_extendedReduceElem_T_2* BD_extendedReduceList_2;
	static BT_extendedReduceElem_T_3* BD_extendedReduceList_3;
	static BT_extendedReduceElem_T* BD_extendedReduceResult_P;
	static BT_extendedReduceElem_T_1* BD_extendedReduceResult_P_1;
	static BT_extendedReduceElem_T_2* BD_extendedReduceResult_P_2;
	static BT_extendedReduceElem_T_3* BD_extendedReduceResult_P_3;
	static PT_bsf_mapElem_T* BD_mapSubList;
	static BT_order_T					BD_order;		// Order for Workers
	static MPI_Status					BD_status[PP_BSF_MAX_MPI_SIZE];		// MPI ststus
	static MPI_Request					BD_request[PP_BSF_MAX_MPI_SIZE];	// MPI request
	static int							BD_sublistSize[PP_BSF_MAX_MPI_SIZE];// SubList Sizes for Workers
	static int							BD_offset[PP_BSF_MAX_MPI_SIZE];		// List offsets for Workers
	/*---------------- Entrance Point ----------------------*/
	void run();
	/*---------------- Assigning Values to BSF Problem Variables ----------------------*/
	void PC_bsfAssignAddressOffset(int value);
	void PC_bsfAssignIterCounter(int value);
	void PC_bsfAssignJobCase(int value);
	void PC_bsfAssignMpiMaster(int value);
	void PC_bsfAssignMpiRank(int value);
	void PC_bsfAssignNumberInSublist(int value);
	void PC_bsfAssignNumOfWorkers(int value);
	void PC_bsfAssignParameter(PT_bsf_parameter_T parameter);
	void PC_bsfAssignSublistLength(int value);
	/*---------------- BSF Skeleton Functions ----------------------*/
	void BC_Init(bool* success);		// Performs the memory allocation and the initialization of the skeleton data structures and variables.
	void BC_Master();				// The head function of the master process.
	void BC_MasterMap(bool exit);	// Forms an order and sends it to the worker processes to perform the Map function in the current iteration.
	void BC_MasterReduce();			// Receives the results from worker processes, collects them in a list, and performs the function Reduce on this list.
	void BC_MpiRun();				// MPI initialization.
	void BC_ProcessExtendedReduceList(BT_extendedReduceElem_T* reduceList, int length, BT_extendedReduceElem_T** extendedReduceResult_P);
	void BC_ProcessExtendedReduceList_1(BT_extendedReduceElem_T_1* reduceList, int length, BT_extendedReduceElem_T_1** extendedReduceResult_P);
	void BC_ProcessExtendedReduceList_2(BT_extendedReduceElem_T_2* reduceList, int length, BT_extendedReduceElem_T_2** extendedReduceResult_P);
	void BC_ProcessExtendedReduceList_3(BT_extendedReduceElem_T_3* reduceList, int length, BT_extendedReduceElem_T_3** extendedReduceResult_P);
	void BC_Worker();				// The head function of a worker process.
	bool BC_WorkerMap();				// Performs the Map function.
	void BC_WorkerReduce();			// Sends to the master process the element that is the sum of all reduce-sublist elements.

	virtual void PC_bsf_CopyParameter(
		PT_bsf_parameter_T parameterIn,
		PT_bsf_parameter_T* parameterOutP
	) = 0;
	virtual void PC_bsf_Init(
		bool* success
	) = 0;
	virtual void PC_bsf_IterOutput(
		PT_bsf_reduceElem_T* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	) = 0;
	virtual void PC_bsf_IterOutput_1(
		PT_bsf_reduceElem_T_1* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	) = 0;
	virtual void PC_bsf_IterOutput_2(
		PT_bsf_reduceElem_T_2* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	) = 0;
	virtual void PC_bsf_IterOutput_3(
		PT_bsf_reduceElem_T_3* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double elapsedTime,
		int nextJob
	) = 0;
	virtual void PC_bsf_JobDispatcher(
		PT_bsf_parameter_T* parameter,
		int* job,
		bool* exit,
		double t
	) = 0;
	//
	virtual void PC_bsf_MapF(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T* reduceElem,
		int* success
	) = 0;
	virtual void PC_bsf_MapF_1(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T_1* reduceElem,
		int* success
	) = 0;
	virtual void PC_bsf_MapF_2(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T_2* reduceElem,
		int* success
	) = 0;
	virtual void PC_bsf_MapF_3(
		PT_bsf_mapElem_T* mapElem,
		PT_bsf_reduceElem_T_3* reduceElem,
		int* success
	) = 0;
	virtual void PC_bsf_ParametersOutput(
		PT_bsf_parameter_T parameter
	) = 0;
	virtual void PC_bsf_ProblemOutput(
		PT_bsf_reduceElem_T* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	) = 0;
	virtual void PC_bsf_ProblemOutput_1(
		PT_bsf_reduceElem_T_1* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	) = 0;
	virtual void PC_bsf_ProblemOutput_2(
		PT_bsf_reduceElem_T_2* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	) = 0;
	virtual void PC_bsf_ProblemOutput_3(
		PT_bsf_reduceElem_T_3* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T parameter,
		double t
	) = 0;
	virtual void PC_bsf_ProcessResults(
		PT_bsf_reduceElem_T* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	) = 0;
	virtual void PC_bsf_ProcessResults_1(
		PT_bsf_reduceElem_T_1* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	) = 0;
	virtual void PC_bsf_ProcessResults_2(
		PT_bsf_reduceElem_T_2* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	) = 0;
	virtual void PC_bsf_ProcessResults_3(
		PT_bsf_reduceElem_T_3* reduceResult,
		int reduceCounter,
		PT_bsf_parameter_T* parameter,
		int* nextJob,
		bool* exit
	) = 0;
	virtual void PC_bsf_ReduceF(
		PT_bsf_reduceElem_T* x,
		PT_bsf_reduceElem_T* y,
		PT_bsf_reduceElem_T* z
	) = 0;
	virtual void PC_bsf_ReduceF_1(
		PT_bsf_reduceElem_T_1* x,
		PT_bsf_reduceElem_T_1* y,
		PT_bsf_reduceElem_T_1* z
	) = 0;
	virtual void PC_bsf_ReduceF_2(
		PT_bsf_reduceElem_T_2* x,
		PT_bsf_reduceElem_T_2* y,
		PT_bsf_reduceElem_T_2* z
	) = 0;
	virtual void PC_bsf_ReduceF_3(
		PT_bsf_reduceElem_T_3* x,
		PT_bsf_reduceElem_T_3* y,
		PT_bsf_reduceElem_T_3* z
	) = 0;
	virtual void PC_bsf_SetListSize(
		int* listSize
	) = 0;
	virtual void PC_bsf_SetInitParameter(
		PT_bsf_parameter_T* parameter
	) = 0;
	virtual void PC_bsf_SetMapListElem(
		PT_bsf_mapElem_T* elem,
		int i
	) = 0;
};


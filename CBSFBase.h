#pragma once
#include "BSF-Types.h"				// Problem Independent Types 
#include "BSF-Data.h"				// Problem Independent Variables & Data Structures 
#include "BSF-SkeletonVariables.h"	// Skeleton Variables

using namespace std;

class CBSFBase
{
public:
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


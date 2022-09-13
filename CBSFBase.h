#pragma once
#include "BSF-Types.h"				// Problem Independent Types 
#include "BSF-Data.h"				// Problem Independent Variables & Data Structures 
#include "BSF-ProblemFunctions.h"	// Predefined Problem Function Forwards
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
#include "BSF-Forwards.h"			// Problem Independent Function Forwards


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
};


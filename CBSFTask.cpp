#include "CBSFTask.h"

void CBSFTask::PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	for (int j = 0; j < PD_n; j++) // Generating initial approximation
		parameter->x[j] = PD_u[j];
};

void CBSFTask::PC_bsf_Init(bool* success) {
	PD_state = PP_STATE_START;
	PD_problemName = PP_PROBLEM_NAME;

	/*debug5*
		if (PP_MODE_BLOCK_HCV_VARIABLE && PP_MODE_USE_LCV_VARIABLE) {
			cout << "Modes PP_MODE_BLOCK_HCV_VARIABLE & PP_MODE_USE_LCV_VARIABLE are incompatible!\n";
			*success = false;
			return;
		}
	/*end debug*/

	*success = LoadMatrixFormat();
	if (*success == false)
		return;

	//
	//
	//

	PD_MTX_File_so = PP_PATH;
	PD_MTX_File_so += PP_MTX_PREFIX;
	PD_MTX_File_so += PD_problemName;
	PD_MTX_File_so += PP_MTX_POSTFIX_SO;

	PD_firstLcvI = INT_MAX;
	PD_firstZcvI = INT_MAX;
	for (int j1 = 0; j1 < PD_n; j1++) {
		if (fabs(PD_c[PD_objI[j1]]) / fabs(PD_c[PD_objI[0]]) <= PP_LOW_COST_PERCENTILE) {
			PD_firstLcvI = j1;
			for (int j2 = j1; j2 < PD_n; j2++)
				if (PD_c[PD_objI[j2]] == 0) {
					PD_firstZcvI = j2;
					PD_minCostPercentile = fabs(PD_c[PD_objI[j2 - 1]]) / fabs(PD_c[PD_objI[0]]);
					break;
				}
			break;
		}
	}

	MakeObjVector(PD_c, PD_objVector);
	ObjUnitVector(PD_unitObjVector);
	Vector_MultiplyByNumber(PD_unitObjVector, PP_DISTANCE_TO_APEX, PD_direction);
	Vector_Addition(PD_u, PD_direction, PD_apexPoint);
}

void CBSFTask::PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_MM;
}

void CBSFTask::PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->a = PD_A[i];
	elem->b = &(PD_b[i]);
}

// 0. Pseudo-pojection
void CBSFTask::PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success // 1 - reduceElem was produced successfully; 0 - otherwise
) {
	int exitCode;

	Vector_ProjectOnHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b, reduceElem->projection, &exitCode);

	switch (exitCode) {
	case PP_EXITCODE_DEGENERATE_INEQUALITY:
		for (int j = 0; j < PD_n; j++)
			reduceElem->projection[j] = 0;
		reduceElem->nonZeroCounter = 0;
		reduceElem->pointIn = true;
		break;
	case PP_EXITCODE_POINT_BELONGS_TO_HALFSPACE:
		for (int j = 0; j < PD_n; j++)
			reduceElem->projection[j] = 0;
		reduceElem->nonZeroCounter = 0;
		reduceElem->pointIn = true;
		break;
	case PP_EXITCODE_SUCCESSFUL_PROJECTING:
		reduceElem->nonZeroCounter = 1;
		reduceElem->pointIn = false;
		break;
	default:
		cout << "Vector_ProjectOnHalfspace: Undefined exit code!" << endl;
		break;
	}
}

// 1. CheckIn
void CBSFTask::PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	if (mapElem->a[PP_ADD_FLAG] == 0)
		reduceElem->pointIn = PointInHalfspace(BSF_sv_parameter.x, mapElem->a, *mapElem->b);
	else
		reduceElem->pointIn = true;
}

void CBSFTask::PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// not used
}

void CBSFTask::PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// not used
}

// 0. Pseudo-pojection
void CBSFTask::PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	Vector_Addition(x->projection, y->projection, z->projection);
	z->nonZeroCounter = x->nonZeroCounter + y->nonZeroCounter;
	z->pointIn = x->pointIn && y->pointIn;
}

// 1. CheckIn
void CBSFTask::PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	z->pointIn = x->pointIn && y->pointIn;
}

void CBSFTask::PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	// not used
}

void CBSFTask::PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// not used
}

//0. Start
void CBSFTask::PC_bsf_ProcessResults(
	PT_bsf_reduceElem_T* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> PC_bsf_ProcessResults: Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = "
			<< PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	};
#endif // PP_MAX_ITER_COUNT

	PD_pointIn = reduceResult->pointIn;
	if (PD_pointIn)
		return;

	static PT_vector_T relaxationVector;
	for (int j = 0; j < PD_n; j++)
		relaxationVector[j] = reduceResult->projection[j] / (double)(reduceResult->nonZeroCounter + 1);
	Vector_PlusEquals(parameter->x, relaxationVector);
}

// 1. Movement on Polytope  ========================================================
void CBSFTask::PC_bsf_ProcessResults_1(
	PT_bsf_reduceElem_T_1* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
#ifdef PP_MAX_ITER_COUNT
	if (BSF_sv_iterCounter > PP_MAX_ITER_COUNT) {
		cout << "-------------> Acceptable maximum number of iterations is exceeded: PP_MAX_ITER_COUNT = " << PP_MAX_ITER_COUNT << endl;
		*exit = true;
		return;
	}
#endif // PP_MAX_ITER_COUNT

	PD_pointIn = reduceResult->pointIn;
}

void CBSFTask::PC_bsf_ProcessResults_2(
	PT_bsf_reduceElem_T_2* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void CBSFTask::PC_bsf_ProcessResults_3(
	PT_bsf_reduceElem_T_3* reduceResult,
	int reduceCounter, // Number of successfully produced Elements of Reduce List
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* nextJob,
	bool* exit // "true" if Stopping Criterion is satisfied, and "false" otherwise
) {
	// not used
}

void CBSFTask::PC_bsf_JobDispatcher(
	PT_bsf_parameter_T* parameter, // Current Approximation
	int* job,
	bool* exit,
	double t
) {
	static PT_vector_T shiftBasePoint;
	static double* ptr_unitVectorToSurface;
	static int landingNo;
	const char* x0_File = PD_MTX_File_x0.c_str();
	//
	bool goOn, repeat;

	static int detDirSwitch;
	static int lcvI;
	static double c_lcvI;
	static double new_c_lcv;
	static int max_lcvI;
	static double max_objF_lcv;
	static double max_new_c_lcv;
	static double objF_lcv;

#define BEGIN_LCV_UTILIZATION	0
#define POSITIVE_LCV			1
#define NEGATIVE_LCV			2
#define POSITIVE_ZCV			3
#define NEGATIVE_ZCV			4
#define END_LCV_UTILIZATION		5

	switch (PD_state) {
	case PP_STATE_START://-------------------------- Start -----------------------------
		assert(!PointInPolytope_s(PD_apexPoint));
		// Preparations for finding feasible point
		Vector_Copy(PD_apexPoint, parameter->x);
		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_FIND_FEASIBLE_POINT;
#ifdef PP_DEBUG
		cout << "--------- Finding feasible point ------------\n";
#endif
		break;
	case PP_STATE_DETERMINE_DIRECTION://------------------------- Determine Direction -----------------------------
		if (!PD_pointIn)
			return;

		if (!PP_MODE_USE_LCV_VARIABLE) {
			/*debug00*/
#ifdef PP_DEBUG
			cout << "w =\t\t";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
				cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
			if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
			cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
#endif
			/*end debug*/
		}

		if (PP_MODE_USE_LCV_VARIABLE) {
			if (PD_firstLcvI == INT_MAX) {
				cout << "Error: The PP_MODE_USE_LCV_VARIABLE=true, but there are no low cost variables!\n";
				*exit = true;
				return;
			}

			switch (detDirSwitch) {
			case BEGIN_LCV_UTILIZATION:
				max_lcvI = 0;
				max_objF_lcv = -INFINITY;

				/*debug00*/
#ifdef PP_DEBUG
				cout << "w =\t\t";
				for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
					cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
				if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
				cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
#endif
				/*end debug*/

				PD_objF_w = ObjF(parameter->x);
				lcvI = PD_firstLcvI;
				c_lcvI = PD_c[PD_objI[lcvI]];
				if (PD_firstLcvI < (PD_firstZcvI == INT_MAX ? PD_n : PD_firstZcvI)) {
					PD_c[PD_objI[lcvI]] = fabs(PD_c[PD_objI[lcvI]]) + fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2);
					new_c_lcv = PD_c[PD_objI[lcvI]];
					MakeObjVector(PD_c, PD_objVector);
					detDirSwitch = POSITIVE_LCV;
				}
				else {
					PD_c[PD_objI[lcvI]] = fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2);
					new_c_lcv = PD_c[PD_objI[lcvI]];
					MakeObjVector(PD_c, PD_objVector);
					detDirSwitch = POSITIVE_ZCV;
				}
				PD_c[PD_objI[lcvI]] = c_lcvI;
				Vector_Copy(PD_u, parameter->x);
				Vector_PlusEquals(parameter->x, PD_objVector);
				return;
			case POSITIVE_LCV:
				/*debug01*
#ifdef PP_DEBUG
				if (lcvI == PD_firstLcvI)
					cout << "---------------- Nonzero low cost variables with '+' ----------------\n";
				cout << "#" << lcvI << "|" << PD_objI[lcvI] << ":\t";
				for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
					cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
				if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
				cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
#endif // PP_DEBUG
				/*end debug*/

				objF_lcv = ObjF(parameter->x);

				if (objF_lcv > PD_objF_w + PP_EPS_ZERO_DIR
					&& objF_lcv > PD_objF_u + PP_EPS_ZERO_DIR
					&& objF_lcv > max_objF_lcv)
				{
					max_objF_lcv = objF_lcv;
					max_lcvI = lcvI;
					max_new_c_lcv = new_c_lcv;
				}

				lcvI++;
				c_lcvI = PD_c[PD_objI[lcvI]];
				if (lcvI < (PD_firstZcvI == INT_MAX ? PD_n : PD_firstZcvI)) {
					PD_c[PD_objI[lcvI]] = fabs(PD_c[PD_objI[lcvI]]) + fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2);
					new_c_lcv = PD_c[PD_objI[lcvI]];
					MakeObjVector(PD_c, PD_objVector);
					PD_c[PD_objI[lcvI]] = c_lcvI;
					Vector_Copy(PD_u, parameter->x);
					Vector_PlusEquals(parameter->x, PD_objVector);
					return;
				}

				lcvI = PD_firstLcvI;
				c_lcvI = PD_c[PD_objI[lcvI]];
				PD_c[PD_objI[lcvI]] = -(fabs(PD_c[PD_objI[lcvI]]) + fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2));
				new_c_lcv = PD_c[PD_objI[lcvI]];
				MakeObjVector(PD_c, PD_objVector);
				PD_c[PD_objI[lcvI]] = c_lcvI;
				Vector_Copy(PD_u, parameter->x);
				Vector_PlusEquals(parameter->x, PD_objVector);
				detDirSwitch = NEGATIVE_LCV;
				return;

			case NEGATIVE_LCV:
				/*debug01*
#ifdef PP_DEBUG
				if (lcvI == PD_firstLcvI)
					cout << "---------------- Nonzero low cost variables with '-' ----------------\n";
				cout << "#" << lcvI << "|" << PD_objI[lcvI] << ":\t";
				for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
					cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
				if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
				cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
#endif // PP_DEBUG
				/*end debug*/

				objF_lcv = ObjF(parameter->x);

				if (objF_lcv > PD_objF_w + PP_EPS_ZERO_DIR
					&& objF_lcv > PD_objF_u + PP_EPS_ZERO_DIR
					&& objF_lcv > max_objF_lcv)
				{
					max_objF_lcv = objF_lcv;
					max_lcvI = lcvI;
					max_new_c_lcv = new_c_lcv;
				}

				lcvI++;
				c_lcvI = PD_c[PD_objI[lcvI]];
				if (lcvI < (PD_firstZcvI == INT_MAX ? PD_n : PD_firstZcvI)) {
					PD_c[PD_objI[lcvI]] = -(fabs(PD_c[PD_objI[lcvI]]) + fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2));
					new_c_lcv = PD_c[PD_objI[lcvI]];
					MakeObjVector(PD_c, PD_objVector);
					PD_c[PD_objI[lcvI]] = c_lcvI;
					Vector_Copy(PD_u, parameter->x);
					Vector_PlusEquals(parameter->x, PD_objVector);
					return;
				}

				if (PD_firstZcvI != INT_MAX) {
					lcvI = PD_firstZcvI;
					c_lcvI = PD_c[PD_objI[lcvI]];
					PD_c[PD_objI[lcvI]] = fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2);
					new_c_lcv = PD_c[PD_objI[lcvI]];
					MakeObjVector(PD_c, PD_objVector);
					PD_c[PD_objI[lcvI]] = c_lcvI;
					Vector_Copy(PD_u, parameter->x);
					Vector_PlusEquals(parameter->x, PD_objVector);
					detDirSwitch = POSITIVE_ZCV;
					return;
				}
				else
				{
					if (max_lcvI > 0) {
						c_lcvI = PD_c[PD_objI[max_lcvI]];
						PD_c[PD_objI[max_lcvI]] = max_new_c_lcv;
						MakeObjVector(PD_c, PD_objVector);
						PD_c[PD_objI[max_lcvI]] = c_lcvI;
					}
					else
						MakeObjVector(PD_c, PD_objVector);

					Vector_Copy(PD_u, parameter->x);
					Vector_PlusEquals(parameter->x, PD_objVector);
					detDirSwitch = END_LCV_UTILIZATION;
					return;
				}
			case POSITIVE_ZCV:
				/*debug01*
#ifdef PP_DEBUG
				if (lcvI == PD_firstZcvI)
					cout << "---------------- Zero cost variables with '+' ----------------\n";
				cout << "#" << lcvI << "|" << PD_objI[lcvI] << ":\t";
				for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
					cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
				if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
				cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
#endif // PP_DEBUG
				/*end debug*/

				objF_lcv = ObjF(parameter->x);

				if (objF_lcv > PD_objF_w + PP_EPS_ZERO_DIR
					&& objF_lcv > PD_objF_u + PP_EPS_ZERO_DIR
					&& objF_lcv > max_objF_lcv)
				{
					max_objF_lcv = objF_lcv;
					max_lcvI = lcvI;
					max_new_c_lcv = new_c_lcv;
				}

				lcvI++;
				if (lcvI < PD_n) {
					c_lcvI = PD_c[PD_objI[lcvI]];
					PD_c[PD_objI[lcvI]] = fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2);
					new_c_lcv = PD_c[PD_objI[lcvI]];
					if (lcvI < PD_firstZcvI)
						PD_c[PD_objI[lcvI]] += c_lcvI;
					MakeObjVector(PD_c, PD_objVector);
					PD_c[PD_objI[lcvI]] = c_lcvI;
					Vector_Copy(PD_u, parameter->x);
					Vector_PlusEquals(parameter->x, PD_objVector);
					return;
				}

				lcvI = PD_firstZcvI;
				c_lcvI = PD_c[PD_objI[lcvI]];
				PD_c[PD_objI[lcvI]] = -fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2);
				new_c_lcv = PD_c[PD_objI[lcvI]];
				MakeObjVector(PD_c, PD_objVector);
				PD_c[PD_objI[lcvI]] = c_lcvI;
				Vector_Copy(PD_u, parameter->x);
				Vector_PlusEquals(parameter->x, PD_objVector);
				detDirSwitch = NEGATIVE_ZCV;
				return;
			case NEGATIVE_ZCV:
				/*debug01*
#ifdef PP_DEBUG
				if (lcvI == PD_firstZcvI)
					cout << "---------------- Zero cost variables with '-' ----------------\n";
				cout << "#" << lcvI << "|" << PD_objI[lcvI] << ":\t";
				for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
					cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
				if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
				cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
#endif // PP_DEBUG
				/*end debug*/

				objF_lcv = ObjF(parameter->x);
				if (objF_lcv > PD_objF_w + PP_EPS_ZERO_DIR
					&& objF_lcv > PD_objF_u + PP_EPS_ZERO_DIR
					&& objF_lcv > max_objF_lcv)
				{
					max_objF_lcv = objF_lcv;
					max_lcvI = lcvI;
					max_new_c_lcv = new_c_lcv;
				}

				lcvI++;
				if (lcvI < PD_n) {
					c_lcvI = PD_c[PD_objI[lcvI]];
					PD_c[PD_objI[lcvI]] = -fabs(PD_c[PD_objI[PD_firstLcvI - 1]] / 2);
					new_c_lcv = PD_c[PD_objI[lcvI]];
					MakeObjVector(PD_c, PD_objVector);
					PD_c[PD_objI[lcvI]] = c_lcvI;
					Vector_Copy(PD_u, parameter->x);
					Vector_PlusEquals(parameter->x, PD_objVector);
					return;
				}

				if (max_lcvI > 0) {
					c_lcvI = PD_c[PD_objI[max_lcvI]];
					PD_c[PD_objI[max_lcvI]] = max_new_c_lcv;
					MakeObjVector(PD_c, PD_objVector);
					PD_c[PD_objI[max_lcvI]] = c_lcvI;
				}
				else
					MakeObjVector(PD_c, PD_objVector);

				Vector_Copy(PD_u, parameter->x);
				Vector_PlusEquals(parameter->x, PD_objVector);
				detDirSwitch = END_LCV_UTILIZATION;
				return;

			case END_LCV_UTILIZATION:
				/*debug02*/
#ifdef PP_DEBUG
				if (max_lcvI > 0) {
					cout << "Optimal LCV is found with native c = "
						<< (PD_c[PD_objI[max_lcvI]] == 0 ? 0 : PD_c[PD_objI[max_lcvI]])
						<< " and new c = " << max_new_c_lcv << " : #" << max_lcvI << "|"
						<< PD_objI[max_lcvI] << ":\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
				}
#endif // PP_DEBUG
				/*end debug*/

				if (max_lcvI == 0)
					cout << "Optimal LCV not found!\n";

				MakeObjVector(PD_c, PD_objVector);
				break;
			default:
				cout << "PC_bsf_JobDispatcher:switch (detDirSwitch): Undefined Switch!" << endl;
				*exit = true;
				return;
			}
		}

		detDirSwitch = BEGIN_LCV_UTILIZATION;

		/*debug00*/
#ifdef PP_DEBUG
		cout << "w =\t\t";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
		if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
		cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x);
		cout << endl;
#endif
		/*end debug*/

		DetermineDirection(parameter, exit, &repeat);

		if (*exit)
			break;

		if (repeat) {
			// Preparations for determining direction
			Vector_Copy(PD_u, parameter->x);
			PD_objF_u = ObjF(PD_u);
			Vector_PlusEquals(parameter->x, PD_objVector);
			PD_numDetDir = 0;
			*job = PP_JOB_PSEUDOPOJECTION;
			PD_state = PP_STATE_DETERMINE_DIRECTION;
#ifdef PP_DEBUG
			cout << "--------- Determine Direction ------------\n";
#endif
			break;
		}

		/*debug*/
		//SavePoint(PD_u, x0_File, t);
		/*end debug*/

		// Preparations for motion
		PD_shiftLength = PP_START_SHIFT_LENGTH;
		PD_numShiftsSameLength = 0;
		Shift(PD_u, PD_direction, PD_shiftLength, parameter->x);
		/*debug00*/
#ifdef PP_DEBUG
		cout << "--------- Moving along surface ------------\n";
		cout << "Sift = " << setw(PP_SETW) << PD_shiftLength << "\tt = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x);
		cout << endl;
#endif 
		/*end debug*/
		* job = PP_JOB_CHECK;
		ptr_unitVectorToSurface = PD_direction;
		PD_state = MOVING_ALONG_SURFACE;
		break;
	case MOVING_ALONG_SURFACE://-------------------------- Moving along surface -----------------------------
		MovingOnSurface(ptr_unitVectorToSurface, PD_u, parameter->x, &goOn);
		if (goOn)
			return;

		// Preparations for landing
		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_LANDING;
		/*debug00*/
#ifdef PP_DEBUG
		cout << "--------- Landing ------------\n";
#endif // 
		/*end debug*/
		break;
	case PP_STATE_LANDING://-------------------------- Landing -----------------------------
		if (!PD_pointIn)
			return;

		Vector_Copy(parameter->x, PD_u);
		PD_objF_u = ObjF(PD_u);

		//WriteTrace(PD_u);

		/*debug8*/
//		if (PD_objF_u + PP_EPS_OBJ > PP_EXACT_OBJ_VALUE) {
		if (fabs(PD_objF_u - PP_EXACT_OBJ_VALUE) <= PP_EPS_OBJ) {
			*exit = true;
			return;
		}
		/*end debug*/


#ifdef PP_DEBUG
		if (landingNo % PP_BSF_TRACE_COUNT == 0) {
			cout << "Landing# " << landingNo << ". Elapsed time: " << round(t) << endl;
			cout << "u =\t\t";
			for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
				cout << setw(PP_SETW) << PD_u[PD_objI[j]];
			if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
			cout << "\tF(t) = " << setw(PP_SETW) << PD_objF_u;
			cout << endl;
			//SavePoint(PD_u, x0_File, t);
		}
		landingNo++;
#endif

		// Preparations for determining direction
		Vector_PlusEquals(parameter->x, PD_objVector);

		/*debug*/if (PointInPolytope_s(parameter->x)) {
			cout << "Apex point in polytope!\n";
			cout << "Maybe, you should increase PP_OBJECTIVE_VECTOR_LENGTH.\n";
			*exit = true;
			return;
		}/*end debug*/

		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
		PD_numDetDir = 0;
		/*debug00*/
#ifdef PP_DEBUG
		cout << "--------- Determining direction ------------\n";
#endif
		/*end debug*/
		break;
	case PP_STATE_FIND_FEASIBLE_POINT://-------------------------- Finding feasible point -----------------------------
		if (!PD_pointIn) {
			return;
		}

		/*debug2*
		SavePoint(parameter->x, x0_File, t);
		cout << "==================> F(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
		*exit = true;
		return;
		/*end debug*/


#ifdef PP_DEBUG
		cout << "Iter # " << BSF_sv_iterCounter << ". Elapsed time: " << round(t) << endl;
		cout << "u =\t\t";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << parameter->x[PD_objI[j]];
		if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";
		cout << "\tF(t) = " << setw(PP_SETW) << ObjF(parameter->x) << endl;
#endif

		Vector_Copy(parameter->x, PD_u);
		PD_objF_u = ObjF(PD_u);

		// Preparations for determining direction
		Vector_PlusEquals(parameter->x, PD_objVector);
		assert(!PointInPolytope_s(parameter->x));

		*job = PP_JOB_PSEUDOPOJECTION;
		PD_state = PP_STATE_DETERMINE_DIRECTION;
		PD_numDetDir = 0;
		/*debug00*/
#ifdef PP_DEBUG
		cout << "--------- Determine Direction ------------\n";
#endif 
		/*end debug*/
		break;
	default://------------------------------------- default -----------------------------------
		cout << "PC_bsf_JobDispatcher: Undefined state!" << endl;
		*exit = true;
		break;
	}
}

void CBSFTask::PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
	cout << "=================================================== BSF Target ====================================================" << endl;
	cout << "Problem name: " << PD_problemName << endl;
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif 
#else
	cout << "OpenMP is turned off!" << endl;
#endif 

#ifdef PP_BSF_FRAGMENTED_MAP_LIST
	cout << "Map List is Fragmented." << endl;
#else
	cout << "Map List is not Fragmented." << endl;
#endif
	cout << "Before conversion: m =\t" << PP_M << "\tn = " << PP_N << endl;
	cout << "After conversion:  m =\t" << PD_m << "\tn = " << PD_n << endl;
	cout << "Eps Zero Compare:\t" << PP_EPS_ZERO_COMPARE << endl;
	cout << "Eps Min Dir Length:\t" << PP_EPS_DIR_LENGTH << endl;
	cout << "Eps Objective:\t\t" << PP_EPS_OBJ << endl;
	cout << "Eps Shift:\t\t" << PP_EPS_SHIFT << endl;
	cout << "Eps Zero Direction:\t" << PP_EPS_ZERO_DIR << endl;
	cout << "Exact Obj Value:\t" << PP_EXACT_OBJ_VALUE << endl;
	cout << "Distance to Apex:\t" << PP_DISTANCE_TO_APEX << endl;
	cout << "Low Cost Percentile:\t" << PP_LOW_COST_PERCENTILE << endl;
	cout << "Gap Max:\t\t" << PP_GAP << endl;
	cout << "Obj Vector Length:\t" << PP_OBJECTIVE_VECTOR_LENGTH << endl;
	cout << "Start Shift Lengt:\t" << PP_START_SHIFT_LENGTH << endl;
	cout << "Blocking obj var:\t" << (PP_MODE_BLOCK_HCV_VARIABLE ? "true" : "false") << endl;
	cout << "Use low cost vars:\t" << (PP_MODE_USE_LCV_VARIABLE ? "true" : "false") << endl;
	cout << "--------------- Statisics ---------------\n";
	cout << "Number of high-cost variables:\t" << (PD_firstLcvI == INT_MAX ? PD_n : PD_firstLcvI) << endl;
	cout << "Number of low-cost variables:\t" << (PD_firstLcvI == INT_MAX ? 0 : (PD_firstZcvI == INT_MAX ? PD_n : PD_firstZcvI) - PD_firstLcvI) << endl;
	cout << "Number of zero-cost variables:\t" << (PD_firstZcvI == INT_MAX ? 0 : PD_n - PD_firstZcvI) << endl;
	cout << "Problem Scale:\t\t\t" << ProblemScale() << endl;
	cout << "--------------- Data ---------------\n";
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")";
		for (int j = 0; j < PD_n; j++)
			cout << setw(PP_SETW) << PD_A[i][j];
		cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
	}
#endif // PP_MATRIX_OUTPUT
	cout << "Obj Function:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_c[PD_objI[j]];
	if (PP_OUTPUT_LIMIT < PD_n)
		cout << " ...";
	cout << endl;
	cout << "Start point:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_u[PD_objI[j]];
	if (PP_OUTPUT_LIMIT < PD_n)
		cout << " ...";
	cout << "\tF(x) = " << ObjF(PD_u);
	cout << endl;
	cout << "Apex point:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_apexPoint[PD_objI[j]];
	if (PP_OUTPUT_LIMIT < PD_n)
		cout << " ...";
	cout << "\tF(x) = " << ObjF(PD_apexPoint);
	cout << endl;
	cout << "-------------------------------------------" << endl;
}

void CBSFTask::PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	for (int j = 0; j < PD_n; j++)
		parameterOutP->x[j] = parameterIn.x[j];

	if (parameterIn.m > 0) {
		PD_A[parameterIn.m][PP_ADD_FLAG] = 1;
		PD_A[parameterIn.m][PD_objI[parameterIn.indexToBlock]] = -parameterIn.sign;
		PD_b[parameterIn.m] = -parameterIn.sign * parameterIn.b;
		PD_A[parameterIn.m + 1][PP_ADD_FLAG] = 1;
		PD_A[parameterIn.m + 1][PD_objI[parameterIn.indexToBlock]] = parameterIn.sign;
		PD_b[parameterIn.m + 1] = parameterIn.sign * parameterIn.b;
	}
}

// 0. Start
void CBSFTask::PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {
	cout << "# " << BSF_sv_iterCounter << "\tTime " << round(elapsedTime);
	cout << "\tx = ";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << parameter.x[PD_objI[j]];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << "\tF(x)= " << setw(PP_SETW) << ObjF(parameter.x) << endl;

	/*debug3*
	const char* x0_File = PD_MTX_File_x0.c_str();
	SavePoint(parameter.x, x0_File, elapsedTime);
	/*end debug*/
}

// 1. Movement on Polytope
void CBSFTask::PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	/*cout << "------------------ 1. Movement on Polytope. Iter # " << BSF_sv_iterCounter << " ------------------" << endl;
	cout << "Elapsed time: " << round(elapsedTime) << endl;
	cout << "PD_u:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_u[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "PD_direction:";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
		cout << setw(PP_SETW) << PD_direction[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "Sift Length = " << PD_shiftLength << endl;/**/
};

// 2.
void CBSFTask::PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}


void CBSFTask::PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// not used
}

// 0. Start
void CBSFTask::PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	ProblemOutput(t);
}

// 1. Movement on Polytope
void CBSFTask::PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	ProblemOutput(t);
}

void CBSFTask::PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

void CBSFTask::PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// not used
}

//---------------------------------- Problem functions -------------------------
inline PT_float_T CBSFTask::Vector_DotProductSquare(PT_vector_T x, PT_vector_T y) {
	PT_float_T sum = 0;
	for (int j = 0; j < PD_n; j++)
		sum += x[j] * y[j];
	return sum;
}

inline PT_float_T CBSFTask::Vector_Norm(PT_vector_T x) {
	return sqrt(Vector_NormSquare(x));
}

inline PT_float_T CBSFTask::Vector_NormSquare(PT_vector_T x) {
	PT_float_T sum = 0;

	for (int j = 0; j < PD_n; j++) {
		sum += x[j] * x[j];
	}
	return sum;
}

inline bool CBSFTask::PointInHalfspace // If the point belongs to the Halfspace with prescigion of PD_Gap
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return Vector_DotProductSquare(a, point) <= b + PP_GAP;
}

inline bool CBSFTask::PointInHalfspace_s // If the point belongs to the Halfspace with prescigion of PP_EPS_ZERO_COMPARE
(PT_vector_T point, PT_vector_T a, PT_float_T b) {
	return Vector_DotProductSquare(a, point) <= b + PP_EPS_ZERO_COMPARE;
}

inline bool CBSFTask::PointInPolytope_s(PT_vector_T x) { // If the point belongs to the polytope with prescigion of PP_EPS_ZERO_COMPARE
	for (int i = 0; i < PD_m; i++) {
		if (PD_A[i][PP_ADD_FLAG] == 1)
			continue;
		if (!PointInHalfspace_s(x, PD_A[i], PD_b[i]))
			return false;
	}
	return true;
}

inline void CBSFTask::Shift(PT_vector_T basePoint, PT_vector_T direction, PT_float_T siftLength, PT_vector_T endPoint) {
	for (int j = 0; j < PD_n; j++)
		endPoint[j] = basePoint[j] + direction[j] * siftLength;
}

// Calculating unit vector of direction from startPoint to endPoint
inline bool CBSFTask::GetDirection(PT_vector_T startPoint, PT_vector_T endPoint, PT_vector_T unitVector) {
	for (int j = 0; j < PD_n; j++) {
		unitVector[j] = endPoint[j] - startPoint[j];
	}
	double normOfUnitVector = Vector_Norm(unitVector);

	if (normOfUnitVector < PP_EPS_DIR_LENGTH)
		return false;

	for (int j = 0; j < PD_n; j++) {
		unitVector[j] /= normOfUnitVector;
	}
	return true;
}

inline void CBSFTask::Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
	for (int j = 0; j < PD_n; j++)
		toPoint[j] = fromPoint[j];
}

inline void CBSFTask::Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
	for (int j = 0; j < PD_n; j++)
		equalVector[j] += plusVector[j];
}

inline void CBSFTask::Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
	for (int j = 0; j < PD_n; j++)
		equalPoint[j] -= minusVector[j];
}

inline void CBSFTask::Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] + y[j];
}

inline void CBSFTask::Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
	for (int j = 0; j < PD_n; j++)
		z[j] = x[j] - y[j];
}

inline void CBSFTask::Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] * r;
}

inline void CBSFTask::Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
	for (int j = 0; j < PD_n; j++)
		x[j] *= r;
}

inline void CBSFTask::Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
	for (int j = 0; j < PD_n; j++)
		x[j] /= r;
}

inline void CBSFTask::Vector_ResetToZero(PT_vector_T x) {  // x = 0
	for (int j = 0; j < PD_n; j++) x[j] = 0;
}

inline void CBSFTask::Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = x/r
	for (int j = 0; j < PD_n; j++)
		y[j] = x[j] / r;
}

inline void CBSFTask::Vector_Round(PT_vector_T x) {
	double floorValue;
	double fractionalPart;
	double sign;
	double absValue;

	for (int j = 0; j < PD_n; j++) {
		if (fabs(x[j]) < PP_EPS_ZERO_DIR) {
			x[j] = 0;
			continue;
		}
		absValue = fabs(x[j]);
		sign = x[j] > 0 ? 1 : -1;
		floorValue = floor(absValue);
		fractionalPart = absValue - floorValue;
		if (1 - fractionalPart < PP_EPS_ZERO_DIR) {
			x[j] = sign * (floorValue + 1);
			continue;
		}
		if (fractionalPart < PP_EPS_ZERO_DIR)
			x[j] = sign * floorValue;
	}
}

inline void CBSFTask::Vector_EpsZero(PT_vector_T x) { // If x[j] < PP_EPS_ZERO_DIR then x[j] = 0
	for (int j = 0; j < PD_n; j++)
		if (fabs(x[j]) < PP_EPS_ZERO_DIR)
			x[j] = 0;
}

inline void CBSFTask::Vector_Unit(PT_vector_T vector) { // Calculating unit vector
	double normOfVector = Vector_Norm(vector);
	for (int j = 0; j < PD_n; j++) {
		vector[j] /= normOfVector;
	}
};

inline PT_float_T CBSFTask::ObjF(PT_vector_T x) {
	PT_float_T s = 0;
	for (int j = 0; j < PD_n; j++)
		s += PD_c[j] * x[j];
	return s;
}

bool CBSFTask::LoadMatrixFormat() {
	int nor,	// Number of matrix rows
		noc,	// Number of matrix columns
		non,	// Number of non-zero elements
		noe;	// Number of equations
	const char* mtxFile;
	FILE* stream;// Input stream
	char str[80] = { '\0' };
	char* chr = str;

	//--------------- Reading A ------------------
	PD_MTX_File_A = PP_PATH;
	PD_MTX_File_A += PP_MTX_PREFIX;
	PD_MTX_File_A += PD_problemName;
	PD_MTX_File_A += PP_MTX_POSTFIX_A;
	mtxFile = PD_MTX_File_A.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d%d", &nor, &noc, &non) < 3) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Unexpected end of file " << mtxFile << endl;
		return false;
	}

	if (nor >= noc) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Number of rows m = " << nor << " must be < " << "Number of columns n = " << noc << "\n";
		return false;
	}

	if (noc != PP_N) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Invalid input data: PP_N must be = " << noc << "\n";
		return false;
	}

	if (nor != PP_M) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Invalid input data: PP_M must be = " << nor << "\n";
		return false;
	}

	PD_m = noe = nor;
	PD_n = noc;

	if (2 * nor + noc > PP_MM) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Invalid input data: number of inequalities m = " << 2 * nor + noc
			<< " must be < PP_MM + 1 =" << PP_MM + 1 << "\n";
		return false;
	}

	for (int k = 0; k < non; k++) {
		int i, j;

		if (fscanf(stream, "%d%d%s", &i, &j, str) < 3) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file'" << mtxFile << "'." << endl;
			return false;
		}

		i -= 1;
		j -= 1;
		if (i < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Negative row index in'" << mtxFile << "'.\n" << endl;
			return false;
		}
		if (j < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Negative column index in'" << mtxFile << "'.\n" << endl;
			return false;
		}
		PD_A[i][j] = strtod(str, &chr);
	}

	/*debug*
for (int i = 0; i < PD_m; i++) {
	for (int j = 0; j < PD_n; j++)
		cout << PD_A[i][j] << " ";
	cout << endl;
}
/*end debug*/

	fclose(stream);

	//--------------- Reading b ------------------
	PD_MTX_File_b = PP_PATH;
	PD_MTX_File_b += PP_MTX_PREFIX;
	PD_MTX_File_b += PD_problemName;
	PD_MTX_File_b += PP_MTX_POSTFIX_B;
	mtxFile = PD_MTX_File_b.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (noe != nor) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int i = 0; i < noe; i++) {
		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_b[i] = strtod(str, &chr);
	}
	fclose(stream);

	/*debug*
for (int i = 0; i < PD_m; i++)
	cout << PD_b[i] << endl;
/*end debug*/

//--------------- Reading lo ------------------
	PD_MTX_File_lo = PP_PATH;
	PD_MTX_File_lo += PP_MTX_PREFIX;
	PD_MTX_File_lo += PD_problemName;
	PD_MTX_File_lo += PP_MTX_POSTFIX_LO;
	mtxFile = PD_MTX_File_lo.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_lo[j] = strtod(str, &chr);
	}

	fclose(stream);

	//--------------- Reading c ------------------
	PD_MTX_File_c = PP_PATH;
	PD_MTX_File_c += PP_MTX_PREFIX;
	PD_MTX_File_c += PD_problemName;
	PD_MTX_File_c += PP_MTX_POSTFIX_C;
	mtxFile = PD_MTX_File_c.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file" << endl;
			return false;
		}
		PD_c[j] = -strtod(str, &chr);
	}
	fclose(stream);

	/*debug*
for (int i = 0; i < PD_n; i++)
	cout << PD_c[i] << " ";
cout << endl;
/*end debug*/

//--------------- Reading hi ------------------
	PD_MTX_File_hi = PP_PATH;
	PD_MTX_File_hi += PP_MTX_PREFIX;
	PD_MTX_File_hi += PD_problemName;
	PD_MTX_File_hi += PP_MTX_POSTFIX_HI;
	mtxFile = PD_MTX_File_hi.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Failure of opening file '" << mtxFile << "'.\n";
		return false;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 1) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file '" << mtxFile << "'." << endl;
			return false;
		}
		PD_hi[j] = strtod(str, &chr);
	}
	fclose(stream);

	bool error = !Conversion();
	if (error) return false;

	SortObjVarI();

	//--------------- Reading x0 ------------------
	PD_MTX_File_x0 = PP_PATH;
	PD_MTX_File_x0 += PP_MTX_PREFIX;
	PD_MTX_File_x0 += PD_problemName;
	PD_MTX_File_x0 += PP_MTX_POSTFIX_X0;
	mtxFile = PD_MTX_File_x0.c_str();
	stream = fopen(mtxFile, "r+b");

	if (stream == NULL) {
		// Generating Coordinates of starting point
		for (int j = 0; j < PD_n; j++)
			PD_u[j] = 0;
		return true;
	}

	SkipComments(stream);
	if (fscanf(stream, "%d%d", &nor, &noc) < 2) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Unexpected end of file'" << mtxFile << "'." << endl;
		return false;
	}
	if (nor != PD_n) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of rows in'" << mtxFile << "'.\n";
		return false;
	}
	if (noc != 1) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout
			<< "Incorrect number of columnws in'" << mtxFile << "'.\n";
		return false;
	}

	for (int j = 0; j < PD_n; j++) {
		if (fscanf(stream, "%s", str) < 0) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout
				<< "Unexpected end of file" << endl;
			return false;
		}
		PD_u[j] = strtod(str, &chr);
	}
	fclose(stream);
	return true;
}

bool CBSFTask::Conversion() { // Transformation to inequalities & dimensionality reduction
	static PT_float_T fvA[PP_MM]; // Free variable coefficients
	static bool Flag[PP_N];		// Flags of free variables to delete
	static int fvEqI;	// Inequality index of free variable
	static bool single;

	for (int jc = 0; jc < PD_n; jc++)
		if (PD_c[jc] == 0 && PD_hi[jc] == PP_INFINITY && PD_lo[jc] == 0)
			Flag[jc] = true;

	for (int jc = 0; jc < PD_n; jc++) {
		if (!Flag[jc])
			continue;
		for (int i = 0; i < PD_m; i++) { // Find corresponding equation
			if (PD_A[i][jc] == 0)
				continue;
			single = true;
			for (int ii = i + 1; ii < PD_m; ii++) // Vertical uniqueness
				if (PD_A[ii][jc] != 0) {
					single = false;
					break;
				}
			if (single)
				fvEqI = i;
			break;
		}

		if (!single) {
			Flag[jc] = false;
		}
		else
			fvA[fvEqI] = PD_A[fvEqI][jc];
	}

	/*debug*
	for (int j = 0; j < PD_n; j++)
		cout << Flag[j] << " ";
	cout << endl;
	/*end debug*/

	/*debug*
	for (int i = 0; i < PD_m; i++)
		cout << fvA[i] << endl;
	/*end debug*/


	static bool PD_delete[PP_MM]; // Rows to delete
	PT_float_T s;

	for (int i = 0; i < PD_m; i++) { // Check inconsistent end degenerate equation
		s = 0;
		for (int j = 0; j < PD_n; j++)
			s += fabs(PD_A[i][j]);
		if (s == 0) {
			if (PD_b[i] != 0) {
				//
				cout
					<< "Inconsistent equation " << i << ": " << s << " = " << PD_b[i] << endl;
				return false;
				PD_delete[i] = true;
			}
		}
	}

	for (int i = 0; i < PD_m; i++) { // Removing degenerate equations
		if (!PD_delete[i]) continue;
		for (int ii = i; ii < PD_m - 1; ii++) {  // Remove null equation
			for (int j = 0; j < PD_n; j++)
				PD_A[ii][j] = PD_A[ii + 1][j];
			fvA[ii] = fvA[ii + 1];
			PD_b[ii] = PD_b[ii + 1];
		}
		PD_m--;
	}

	for (int jc = 0; jc < PD_n; jc++) { // Delete free variables
		if (!Flag[jc])
			continue;
		for (int j = jc; j < PD_n - 1; j++) { // Delete column
			PD_c[j] = PD_c[j + 1];
			PD_lo[j] = PD_lo[j + 1];
			PD_hi[j] = PD_hi[j + 1];
			Flag[j] = Flag[j + 1];
			for (int i = 0; i < PD_m; i++)
				PD_A[i][j] = PD_A[i][j + 1];
		}

		PD_n--;
		jc--;
		for (int i = 0; i < PD_m; i++)
			PD_A[i][PD_n] = 0;
		PD_c[PD_n] = 0;
		PD_lo[PD_n] = 0;
		PD_hi[PD_n] = 0;
	}

	/*debug*
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")\t";
		for (int j = 0; j < PD_n; j++)
			cout << PD_A[i][j] << " ";
		cout << endl;
	}
	cout << "----------------------------------------\n";
	/*end debug*/

	int m = PD_m;
	for (int i = 0; i < m; i++) { // Conversion to inequalities

		if (fvA[i] == 0) { // Equation without free variable => adding inequality.
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = -PD_A[i][j];
			PD_b[PD_m] = -PD_b[i];
			PD_m++;
			assert(PD_m <= PP_MM);
		}
		else {
			if (fvA[i] < 0) { // Free variable is negative => change sign to opposite.
				for (int j = 0; j < PD_n; j++)
					PD_A[i][j] = -PD_A[i][j];
				PD_b[i] = -PD_b[i];
			}
		}
	}

	/*debug*
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")\t";
		for (int j = 0; j < PD_n; j++)
			cout << PD_A[i][j] << " ";
		cout << endl;
	}
	cout << "----------------------------------------\n";
	/*end debug*/

	for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
		for (int j = 0; j < PD_n; j++)
			if (PD_A[i][j] == 0)
				PD_A[i][j] = 0;

	for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
		for (int j = 0; j < PD_n; j++)
			PD_A[i + PD_m][j] = 0;
		PD_A[i + PD_m][i] = -1;
		if (PD_lo[i] == 0)
			PD_b[i + PD_m] = 0;
		else
			PD_b[i + PD_m] = -PD_lo[i];
	}
	PD_m += PD_n;
	assert(PD_m <= PP_MM);

	for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
		if (PD_hi[i] != PP_INFINITY) {
			for (int j = 0; j < PD_n; j++)
				PD_A[PD_m][j] = 0;
			PD_A[PD_m][i] = 1;
			PD_b[PD_m] = PD_hi[i];
			PD_m++;
			assert(PD_m <= PP_MM);
		}
	}

	/*debug*
	for (int i = 0; i < PD_m; i++) {
		cout << i << ")\t";
		for (int j = 0; j < PD_n; j++)
			cout << PD_A[i][j] << " ";
		cout << endl;
	}
	cout << "----------------------------------------\n";
	/*end debug*/

	/*debug*
	if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
		for (int j = 0; j < PD_n; j++)
			cout << PD_c[j] << endl;
	/*end debug*/

	return true;
}

bool CBSFTask::SavePoint(PT_vector_T x, const char* filename, double elapsedTime) {
	FILE* stream;

	stream = fopen(filename, "w");
	if (stream == NULL) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Failure of opening file '" << filename << "'.\n";
		return false;
	}

	fprintf(stream, "%c Elapsed time: %.0f\n%d %d\n", '%', elapsedTime, PD_n, 1);

	for (int j = 0; j < PD_n; j++)
		fprintf(stream, "%.14f\n", x[j]);

	fclose(stream);
	return true;
}

// Point projection onto Half-space <a,x> <= b
inline void CBSFTask::Vector_ProjectOnHalfspace(PT_vector_T point, PT_vector_T a, PT_float_T b, PT_vector_T projection, int* exitCode) {
	double factor;
	double aNormSquare = Vector_NormSquare(a);

	if (sqrt(aNormSquare) < PP_EPS_ZERO_COMPARE) {
		*exitCode = PP_EXITCODE_DEGENERATE_INEQUALITY;
		return;
	}

	factor = (b - Vector_DotProductSquare(point, a)) / aNormSquare;

	if (factor > 0 || fabs(factor) < PP_EPS_ZERO_COMPARE) {
		*exitCode = PP_EXITCODE_POINT_BELONGS_TO_HALFSPACE;
		return;
	}

	for (int j = 0; j < PD_n; j++)
		projection[j] = factor * a[j];

	*exitCode = PP_EXITCODE_SUCCESSFUL_PROJECTING;
	return;
}

inline PT_float_T CBSFTask::Distance(PT_vector_T x, PT_vector_T y) {
	PT_vector_T z;
	Vector_Subtraction(x, y, z);
	return Vector_Norm(z);
}

inline void CBSFTask::ObjUnitVector(PT_vector_T objUnitVector) { // Calculating Objective Unit Vector
	double c_norm = Vector_Norm(PD_c);
	Vector_DivideByNumber(PD_c, c_norm, objUnitVector);
}

inline void CBSFTask::ShrinkUnitVector(PT_vector_T objUnitVector, int shrinkBound) { // Shrink Objective Unit Vector from 0 to (shrinkBound-1)
	for (int j = 0; j < shrinkBound; j++)
		objUnitVector[PD_objI[j]] = 0;
	for (int j = shrinkBound; j < PD_n; j++)
		objUnitVector[j] = PD_c[PD_objI[j]];
	double norm = Vector_Norm(objUnitVector);
	Vector_DivideEquals(objUnitVector, norm);
}

inline void CBSFTask::MakeObjVector(PT_vector_T c, PT_vector_T objVector) { // Calculating Objective Vector
	double c_norm = Vector_Norm(c);
	Vector_MultiplyByNumber(c, PP_OBJECTIVE_VECTOR_LENGTH / c_norm, objVector);
}

inline void CBSFTask::ProblemOutput(double elapsedTime) {
	cout << "=============================================" << endl;
	cout << "Elapsed time: " << elapsedTime << endl;
	cout << "Iterations: " << BSF_sv_iterCounter << endl;
	cout << "Optimal objective value: " << setw(PP_SETW) << PD_objF_u << endl;
	cout << "Exact objective value:   " << setw(PP_SETW) << PP_EXACT_OBJ_VALUE << endl;
	cout << "Relative error = " << fabs(PD_objF_u - PP_EXACT_OBJ_VALUE) / fabs(PP_EXACT_OBJ_VALUE) << endl;
	cout << "=============================================" << endl;
	//
	const char* solutionFile = PD_MTX_File_so.c_str();
	Vector_EpsZero(PD_u);
	if (SavePoint(PD_u, solutionFile, elapsedTime))
		cout << "Solution is saved into the file '" << solutionFile << "'." << endl;
	cout << "Solution:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_u[j];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
	cout << "Ordered:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_u[PD_objI[j]];
	if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	cout << endl;
}

inline void CBSFTask::SkipComments(FILE* stream) {
	fpos_t pos;	// Position in the input stream
	int res;
	res = fscanf(stream, "\n");
	fgetpos(stream, &pos);
	while (getc(stream) == '%') {
		while (getc(stream) != 10);
		res = fscanf(stream, "\n");
		fgetpos(stream, &pos);
	};
	fsetpos(stream, &pos);
};

inline void CBSFTask::SortObjVarI() { // Sorting objective variables in absolute descending order
	PT_float_T bigestAbsVal;
	int iBig;

	for (int j = 0; j < PD_n; j++)
		PD_objI[j] = j;

	for (int J = 0; J < PD_n - 1; J++) {
		bigestAbsVal = fabs(PD_c[PD_objI[J]]);
		iBig = J;
		for (int j = J + 1; j < PD_n; j++) {
			if (fabs(PD_c[PD_objI[j]]) > bigestAbsVal) {
				bigestAbsVal = fabs(PD_c[PD_objI[j]]);
				iBig = j;
			}
		}
		if (iBig != J) {
			int j;
			j = PD_objI[iBig];
			PD_objI[iBig] = PD_objI[J];
			PD_objI[J] = j;
		}
	}
}

inline void CBSFTask::MovingOnSurface(PT_vector_T ptr_unitVectorToSurface, PT_vector_T basePoint, PT_vector_T x, bool* goOn)
{
	static PT_float_T objF_basePoint = ObjF(basePoint);
	PT_float_T objF_x = ObjF(x);

	//if (PD_pointIn && objF_basePoint < objF_x) {
	if (PD_pointIn) {
		PD_numDetDir = 0;
		PD_numShiftsSameLength++;
		if (PD_numShiftsSameLength > PP_MAX_NUM_SHIFTS_SAME_LENGTH) {
			PD_shiftLength *= 2;
			PD_numShiftsSameLength = 0;
		}
		Vector_Copy(x, basePoint);
		//PD_oneStepDone = true;
		objF_basePoint = objF_x;
		Shift(basePoint, ptr_unitVectorToSurface, PD_shiftLength, x);

		/*debug00*/
#ifdef PP_DEBUG
		cout << "Sift = " << setw(PP_SETW) << PD_shiftLength << "\tt = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << x[PD_objI[j]];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << "\tF(t) = " << setw(PP_SETW) << ObjF(x);
		cout << endl;
#endif 
		/*end debug*/

		* goOn = true;
		return;
	}

	if (PD_shiftLength < PP_EPS_SHIFT) {
		*goOn = false;
		return;
	}

	PD_shiftLength /= 2;
	PD_numShiftsSameLength = 0;
	Shift(basePoint, ptr_unitVectorToSurface, PD_shiftLength, x);
	*goOn = true;
	/**#ifdef PP_DEBUG
		cout << "Sift = " << setw(PP_SETW) << PD_shiftLength << "\tt = ";
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)
			cout << setw(PP_SETW) << x[PD_objI[j]];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
		cout << "\tF(t) = " << setw(PP_SETW) << ObjF(x);
		cout << endl;
	#endif // PP_DEBUG /**/
}

inline void CBSFTask::DetermineDirection(PT_bsf_parameter_T* parameter, bool* exit, bool* repeat) {

	*repeat = false;

	PD_numDetDir++;
	if (PD_numDetDir > 2) {
		cout << "Number of Squential States 'Determine Direction' is greater than 2!\n";
		cout << "Maybe, you should increase PP_OBJECTIVE_VECTOR_LENGTH.\n";
		*exit = true;
		return;
	}

	if (!GetDirection(PD_u, parameter->x, PD_direction)) {
		cout << "Direction is too small!\n";
		cout << "Maybe, you should increas PP_OBJECTIVE_VECTOR_LENGTH.\n";
		*exit = true;
		return;
	}

	/*debug7*
if (PD_objF_u >= ObjF(parameter->x) + PP_EPS_OBJ) {
	cout << setw(PP_SETW) << "F(u) = " << PD_objF_u << " >= F(w) = " << ObjF(parameter->x) << endl;
	if (PP_MODE_BLOCK_HCV_VARIABLE)
		cout << "Maybe, you should make #define PP_MODE_BLOCK_HCV_VARIABLE false.\n";
	else
		cout << "Maybe, you should decreas PP_EPS_OBJ.\n";
	*exit = true;
	return;
}
	/*end debug*/

	/*debug8*
if (fabs(ObjF(parameter->x) - PD_objF_u) < PP_EPS_OBJ) {
	cout << setw(PP_SETW) << "F(u) = " << PD_objF_u << " == F(w) = " << ObjF(parameter->x) << "\n";
	cout << "Maybe, you should decreas PP_EPS_OBJ.\n";
	*exit = true;
	return;
}
	/*end debug*/

	Vector_EpsZero(PD_direction);

	/*debug00*/
#ifdef PP_DEBUG //----------------------------------------------//
	cout << "D =\t\t";											//
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++)		//
		cout << setw(PP_SETW) << PD_direction[PD_objI[j]];		//
	if (PP_OUTPUT_LIMIT < PD_n) cout << " ...";					//
	cout << endl;												//
#endif // PP_DEBUG ---------------------------------------------//
/*end debug*/

	if (PP_MODE_BLOCK_HCV_VARIABLE) {
		parameter->sign = (PD_c[PD_objI[PD_indexToBlock]] >= 0 ? 1 : -1);
		if (((parameter->sign * PD_direction[PD_objI[PD_indexToBlock]] <= 0)
			|| (fabs(PD_direction[PD_objI[PD_indexToBlock]]) < PP_EPS_ZERO_DIR))
			&& PD_indexToBlock < PD_n - 1)
		{
			parameter->indexToBlock = PD_indexToBlock;
			parameter->m = PD_m;
			parameter->b = PD_u[PD_objI[PD_indexToBlock]];
			PD_m += 2;
			assert(PD_m <= PP_MM);
			PD_numDetDir = 0;
			PD_indexToBlock++;
			ShrinkUnitVector(PD_unitObjVector, PD_indexToBlock);
			Vector_MultiplyByNumber(PD_unitObjVector, PP_OBJECTIVE_VECTOR_LENGTH, PD_objVector);
#ifdef PP_DEBUG
			cout << "Variable " << PD_indexToBlock - 1 << " is blocked.\n";
#endif
			* repeat = true;
		}
	}
}

inline double CBSFTask::ProblemScale() {
	double problemScale = 0;
	for (int i = 0; i < PD_m; i++) {
		for (int j = 0; j < PD_n; j++) {
			if (fabs(PD_A[i][j]) < PP_EPS_ZERO_COMPARE)
				continue;
			problemScale = PF_MAX(problemScale, fabs(PD_b[i] / PD_A[i][j]));
		}
	}
	return problemScale;
}

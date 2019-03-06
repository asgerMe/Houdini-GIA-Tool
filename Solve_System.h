#include "GETFEMMatrices.h"
#include "SIM_FEMSolver.h"

#ifndef __SOLVE_SYSTEM__
#define __SOLVE_SYSTEM__


class SOLVE_LINEAR
{
public:
	static SIM_SingleSolver::SIM_Result SOLVE_3D(SIM_Data_Sparse *sparse, GU_Detail *gdp, Material &parms, Wilson &w, unsigned int sub_steps, int solve_explicit);
	static void invert_sparse(GU_Detail *gdp, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &w);
};



#endif
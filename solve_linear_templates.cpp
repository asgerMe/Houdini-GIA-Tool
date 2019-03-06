#include "Solve_System.h"
#include "GetFEMMatrices.h"

template void SOLVE_LINEAR::invert_sparse<Data_struct>(GU_Detail *gdp, Data_struct &ds, const Wilson &w);
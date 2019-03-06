#include "Solve_System.h"
#include "Preconditioner.h"
#include "IterationTester.h"
#include <UT/UT_SparseMatrix.h>
#include <UT/UT_MatrixSolver.h>
#include <SIM/SIM_Engine.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_SingleSolver.h>
#include <SIM/SIM_GeometryCopy.h>
#include <iostream>
#include "SIM_FEMSolver.h"

SIM_SingleSolver::SIM_Result SOLVE_LINEAR::SOLVE_3D(SIM_Data_Sparse *sparse, GU_Detail *gdp, Material &parms, Wilson &w, unsigned int sub_steps, int __explicit__)
{
	Data_struct ds;
	GA_RWHandleV3D wi(gdp->findAttribute(GA_ATTRIB_POINT, "wt"));
	GA_RWHandleD gflux(gdp->findAttribute(GA_ATTRIB_POINT, "G_flux"));
	GA_RWHandleV3 grav(gdp->findAttribute(GA_ATTRIB_POINT, "gravity_force"));
	GA_RWHandleV3D int_n(gdp->findAttribute(GA_ATTRIB_POINT, "interface_normal"));
	GA_RWHandleD field_divergence(gdp->findAttribute(GA_ATTRIB_POINT, "field_divergence"));

	unsigned int n_points = gdp->getNumPoints();
	if (__explicit__ != -2)
	{
		ds.init(n_points, false);
		if (!sparse->is_loaded())
			sparse->setDivisions(3*n_points);
	}
	else
	{
		ds.init(n_points, true);
		if (!sparse->is_loaded())
			sparse->setDivisions(n_points);
	}
	
	if (sparse->is_loaded())
		ds.stashed = true;
	else
		ds.stashed = false;

	if (ds.stashed)
	{
		ds.fixed_global_force = *sparse->fixed_global;
	}

	UT_AutoInterrupt progress0("Substepping");

	for (unsigned int i = 0; i < sub_steps; i++)
	{
		if (gflux.isValid())
			gflux.clear();
		if (wi.isValid())
			wi.clear();
		if (grav.isValid())
			grav.clear();
		if (int_n.isValid())
			int_n.clear();
		if (field_divergence.isValid())
			field_divergence.clear();

		GA_Iterator primit(gdp->getPrimitiveRange());
		GA_PrimitiveGroup *active_prims = gdp->findPrimitiveGroup("FEM_PRIMS");
		ds.i = i;
		bool use_vic = false;
		UT_AutoInterrupt progress("Assembling linear system of equations");
		UT_Vector3D w0(0.0);
		if (progress.wasInterrupted())
			break;
		std::cout << "works" << std::endl;
		for (primit; !primit.atEnd(); ++primit)
		{
			if (progress.wasInterrupted())
				break;
			GA_Offset primoff = *primit;
			const GA_Primitive *prim = gdp->getPrimitive(primoff);
			if (prim->getPointRange().getEntries() != 4)
				return SIM_SingleSolver::SIM_SOLVER_FAIL;

			bool attribute_present = false;

			if (__explicit__ != -2)
			{
				gridUpdate::update_material_parms(gdp, prim, ds, parms, attribute_present, use_vic);
				gridUpdate::assemble_element_matrices(gdp, prim, ds, parms, __explicit__);
				gridUpdate::update_jacobians(gdp, primoff, ds);
				gridUpdate::assemble_B(ds);
				gridUpdate::global_indices(gdp, ds);
			}
			
			if (__explicit__ == -2)  //Solve laplacian 
			{
				if (active_prims == 0 || !ds.stashed)
				{
					gridUpdate::assemble_element_matrices_laplacian(gdp, prim, ds, parms);
					gridUpdate::update_jacobians(gdp, primoff, ds);
					gridUpdate::assemble_B(ds);
					gridUpdate::assemble_B_laplacian(ds);
					gridUpdate::assemble_global_laplacian(gdp, primoff, ds, sparse, w);
				}
				if (active_prims != 0 && ds.stashed)
				{
					if (active_prims->contains(prim))
					{
						gridUpdate::assemble_element_matrices_laplacian(gdp, prim, ds, parms);
						gridUpdate::update_jacobians(gdp, primoff, ds);
						gridUpdate::assemble_B(ds);
						gridUpdate::assemble_B_laplacian(ds);
						gridUpdate::assemble_global_laplacian(gdp, primoff, ds, sparse, w);
					}
					else
						continue;
				}
			}
			
			if (__explicit__ == 0)  //Assemble for implicit update
			{
				gridUpdate::assemble_global_implicit(gdp, primoff, ds, sparse, parms, w, use_vic);
			}
			if (__explicit__ == 1)
			{
				gridUpdate::assemble_global(gdp, primoff, ds, sparse, parms, use_vic);
			}
		}
 		
		if (__explicit__ == -2)
		{
			if (!ds.stashed)
				gridUpdate::set_diagonal(ds, sparse,  gdp);
			ds.global_force.addScaledVec(1, ds.fixed_global_force);
			invert_sparse(gdp, ds, sparse, w);
			gridUpdate::laplacian_solution(gdp, ds);
		}
		
		if (__explicit__ == 1)
		{
			gridUpdate::step_solution_exp(gdp, ds, parms, w);
		}

		if (__explicit__ == 0)
		{
			
			gridUpdate::apply_inertia(ds, sparse, gdp, w);
			invert_sparse(gdp, ds, sparse, w);
			gridUpdate::step_solution_imp(gdp, ds, w, parms);

		}
	}

	if (!ds.stashed)
	{
		*sparse->fixed_global = ds.fixed_global_force;
	}
	return SIM_SingleSolver::SIM_SOLVER_SUCCESS;
};

void SOLVE_LINEAR::invert_sparse(GU_Detail *gdp, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &w)
{
	//solve systems of equations
	
	DiagonalPreconditioner diagonal_preconditioner( *sparse->stiffness_matrix);
	UT_Functor2<void, const UT_Vector&, UT_Vector&> multiply_full(
		&(*sparse->stiffness_matrix), &UT_SparseMatrix::multVec
	);

	int cp = 200;
	double error = 0.000000000000001;
	IterationTester iteration_tester(cp, error);
	UT_Functor2<bool, int, const UT_Vector&> keep_iterating(iteration_tester);
	UT_Functor2<void, const UT_Vector&, UT_Vector&>
		preconditioner_full(diagonal_preconditioner);

	UT_MatrixIterSolver::PCG(
		ds.X, ds.global_force,
		multiply_full,
		preconditioner_full,
		keep_iterating
	);
}






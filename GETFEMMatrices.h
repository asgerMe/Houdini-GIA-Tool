
#include <UT/UT_Vector.h>
#include <UT/UT_Matrix.h>
#include <GU/GU_Detail.h>
#include "Placeholder.h"
#include "MaterialProperties.h"
#include "Collision.h"
#include "SIM_FEMSolver.h"


#ifndef ASSEMBLE_M2H
#define ASSEMBLE_M2H

class gridUpdate
{
public:
	gridUpdate();
	static void update_tetrapoints(GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds, const unsigned int i);

	static void update_jacobians(GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds);

	static void assemble_B(Data_struct &ds);
	static void assemble_B_laplacian(Data_struct &ds);

	//static void assemble_BNL(Data_struct &ds);

	static void assemble_global( GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds, SIM_Data_Sparse *sparse, const Material &parms, const bool &use_viscosity);
	static void assemble_global_laplacian(GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &wil);
	static void assemble_global_implicit(GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds, SIM_Data_Sparse *sparse, const Material &parms, const Wilson &w, const bool &use_viscosity);
	static void inject_element_matrix_laplacian(const GU_Detail *gdp, unsigned int i, unsigned int j, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &wil);
	static void set_diagonal(Data_struct &ds, SIM_Data_Sparse *sparse, const GU_Detail *gdp);
	static void get_surface_force(const GU_Detail *gdp, const GA_Offset &primoff, const UT_Vector4i &tetra_points, Data_struct &ds);
	static void laplacian_solution(GU_Detail *gdp, Data_struct &ds);
    static void assemble_m1(Data_struct &ds, const Material &parms);
	static void assemble_m2(Data_struct &ds, const Material &parms);
	static void assemble_E(Data_struct &ds, const Material &parms);
	static void add_full_stress_vector(Data_struct &ds, GU_Detail *gdp, const GA_Primitive *prim, const bool &use_vic);
	static void update_material_parms(const GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds, Material &parms, bool &attribute_present, bool &use_vic);
	static void form_stress_tensor(Data_struct &ds);
	static void assemble_linear_element_matrix(Data_struct &ds, bool use_vic);
	static void inject_element_matrix(const GU_Detail *gdp, unsigned int i, unsigned int j, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &wil);
	static void step_solution_imp(GU_Detail *gdp, Data_struct &ds, const Wilson &wil, const Material &material);
	static void step_solution_exp(GU_Detail *gdp, Data_struct &ds, Material &parms, const Wilson &wil);

	static double point_to_prim(const GA_Primitive *prim, GA_ROHandleD &attrib_handle, const Data_struct &ds);
	static void assemble_element_matrices(GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds, const Material &parms, int __explicit__);
	static void assemble_element_matrices_laplacian(GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds, const Material &parms);
	static void rotate_degree_of_freedom(GU_Detail *gdp, GA_Offset primoff, Data_struct &ds, const unsigned int &i, const bool stashed);
	static void rotate_degree_of_freedom(const GA_Detail *gdp, Data_struct &ds, const GA_Offset &ptoff, const GA_Index &idx, UT_VectorD &target_vector, bool transpose);

	static void global_indices(GU_Detail *gdp, Data_struct &ds);
	static void apply_inertia(Data_struct &ds, SIM_Data_Sparse *sparse, GU_Detail *gdp, const Wilson &wil);
	static void calculate_and_apply_winkler(const GA_Detail *gdp, double &surface_force, const Data_struct &ds, const unsigned int &i);
	static void calculate_and_apply_winkler(const GA_Detail *gdp, UT_Vector3D &surface_force, const Data_struct &ds, const unsigned int &i);
	
	static void add_winkler(const GA_Detail *gdp, double c, const UT_Vector3D &normal, Data_struct &ds, const unsigned int &i);

};

#endif // !ASSEMBLE_M2H
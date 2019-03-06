#include <UT/UT_Matrix.h>
#include <UT/UT_SparseMatrix.h>
#include <UT/UT_Vector.h>
#include <GU/GU_Detail.h>



#ifndef __PLACEHOLDER__
#define __PLACEHOLDER__

class Data_struct 
{
public:
	void init(unsigned int points, bool laplacian);
	
	double determinant;
	unsigned int i;
	bool stashed;

	UT_MatrixD local_stiffness;
	UT_Matrix4D jacobian;
	UT_Matrix4D inverse_jacobian;
	
	UT_MatrixD B;
	UT_MatrixD B_temp;
	UT_MatrixD BT;

	UT_MatrixD Bs;
	UT_MatrixD BTs;

	UT_Vector3D interface_normal;
	UT_MatrixD pre_stress;
	UT_MatrixD element_stiffness_matrix_nl;
	UT_MatrixD element_stiffness_matrix;
	UT_MatrixD normal_matrix;

	UT_MatrixD non_linear_element_3D;
	mutable UT_MatrixD N;
	UT_MatrixD stress_strain_matrix;

	UT_MatrixD sigma;

	UT_MatrixD M1;
	UT_MatrixD M2;
	UT_VectorD global_displacements;
	UT_VectorD m1_holder;
	UT_VectorD m2_holder;
	UT_VectorD point_values;
	UT_Matrix3D stress;
	UT_VectorD stress_vector;
	UT_VectorD strain_vector;

	UT_Vector4i tetra_points;
	UT_Vector4i tetra_points_offset;
	UT_Vector4i bound;
	UT_Vector4i winkler;
	UT_Vector4i core;
	UT_VectorD vector_displacements;
	fpreal surf_area;
	UT_SparseMatrixD global_stiffness;
	UT_VectorD global_force;
	UT_VectorD fixed_global_force;
	UT_VectorD global_index;
	UT_VectorD global_bound_index;
	UT_VectorD divergence_vector;

	UT_VectorD lumped_mass;
	UT_VectorD X;
};

#endif
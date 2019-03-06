#include "Placeholder.h"
#include <GU/GU_Detail.h>
#include <UT/UT_SparseMatrix.h>
#include <iostream>



void Data_struct::init(unsigned int points, bool laplacian)
{
	determinant = 0;
	stashed = false;
	tetra_points(0) = -1;
	tetra_points(1) = -1;
	tetra_points(2) = -1;
	tetra_points(3) = -1;

	tetra_points_offset(0) = -1;
	tetra_points_offset(1) = -1;
	tetra_points_offset(2) = -1;
	tetra_points_offset(3) = -1;

	bound(0) = 0;
	bound(1) = 0;
	bound(2) = 0;
	bound(3) = 0;

	winkler(0) = 0;
	winkler(1) = 0;
	winkler(2) = 0; 
	winkler(3) = 0;

	core(0) = 0;
	core(1) = 0;
	core(2) = 0;
	core(3) = 0;

	i = 0;
	surf_area = 0;
	interface_normal *=0;


	B.init(0, 5, 0, 11);
	B.zero();

	BT.init(0, 11, 0, 5);
	BT.zero();

	N.init(0, 5, 0, 11);
	N.zero();

	N(0, 0) = 0.25;
	N(0, 3) = 0.25;
	N(0, 6) = 0.25;
	N(0, 9) = 0.25;

	N(1, 1) = 0.25;
	N(1, 4) = 0.25;
	N(1, 7) = 0.25;
	N(1, 10) = 0.25;

	N(2, 2) = 0.25;
	N(2, 5) = 0.25;
	N(2, 8) = 0.25;
	N(2, 11) = 0.25;

	point_values.init(0, 11);
	point_values.zero();

	vector_displacements.init(0, 11);
	vector_displacements.zero();

	divergence_vector.init(0, 11);
	divergence_vector.zero();

	pre_stress.init(0, 11, 0, 11);
	pre_stress.zero();

	if (!laplacian)
	{
		global_force.init(0, 3 * (points)-1);
		global_force.zero();

		fixed_global_force.init(0, 3 * (points)-1);
		fixed_global_force.zero();

		X.init(0, 3 * (points)-1);
		X.zero();

		global_stiffness.init(3 * (points), 3 * (points));
		global_stiffness.zero();

		global_displacements.init(0, 3* (points)-1);
		global_displacements.zero();
		
		sigma.init(0, 5, 0, 5);
		sigma.zero();

		M1.init(0, 5, 0, 5);
		M1.zero();

		M2.init(0, 5, 0, 5);
		M2.zero();

		m1_holder.init(0, 5);
		m1_holder.zero();

		m2_holder.init(0, 5);
		m2_holder.zero();

		element_stiffness_matrix_nl.init(0, 11, 0, 11);
		element_stiffness_matrix_nl.zero();

		element_stiffness_matrix.init(0, 11, 0, 11);
		element_stiffness_matrix.zero();

		normal_matrix.init(0, 11, 0, 11);
		normal_matrix.zero();

		non_linear_element_3D.init(0, 11, 0, 11);
		non_linear_element_3D.zero();

		lumped_mass.init(0, points - 1);
		lumped_mass.zero();

		global_index.init(0, 11);
		global_index.zero();

		stress_vector.init(0, 5);
		stress_vector.zero();

		global_bound_index.init(0, 11);
		global_bound_index.zero();

		stress_strain_matrix.init(0, 8, 0, 8);
		stress_strain_matrix.zero();

	}
	if (laplacian)
	{

		global_force.init(0, (points)-1);
		global_force.zero();

		fixed_global_force.init(0, (points)-1);
		fixed_global_force.zero();


		X.init(0, (points)-1);
		X.zero();

		global_stiffness.init((points), (points));
		global_stiffness.zero();

		Bs.init(0, 4, 0, 4);
		Bs.zero();

		BTs.init(0, 4, 0, 2);
		BTs.zero();

		non_linear_element_3D.init(0, 4, 0, 4);
		non_linear_element_3D.zero();

	}
}





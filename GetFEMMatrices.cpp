#include "GETFEMMatrices.h"
#include "Placeholder.h"
#include "MaterialProperties.h"
#include "Collision.h"
#include "SIM_FEMSolver.h"
#include <SIM/SIM_Object.h>
#include <SIM/SIM_ObjectArray.h>



gridUpdate::gridUpdate()
{

}

void gridUpdate::update_tetrapoints(GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds , const unsigned int i)
{
	GA_RWHandleV3D sf(gdp->findAttribute(GA_ATTRIB_POINT, "surface_force"));
	GA_RWHandleD smd(gdp->findAttribute(GA_ATTRIB_POINT, "surface_mass_density"));
	GA_RWHandleI winklerp(gdp->findAttribute(GA_ATTRIB_POINT, "winkler"));
	GA_ROHandleI corep(gdp->findAttribute(GA_ATTRIB_POINT, "core"));
	UT_Vector3 zero(0.0);

	if (i < 4)
	{
		GA_Index index = prim->getPointIndex(i);
		ds.tetra_points(i) = index;
		UT_Vector4i offsets;
		GA_Offset ptoff = gdp->getPointMap().offsetFromIndex(GA_Index(ds.tetra_points(i)));
		ds.tetra_points_offset(i) = ptoff;
		bool is_surface = Collider::is_surface(gdp, ptoff, "surface_points");

		if (is_surface)
			ds.bound(i) = 1;
		else
		{
			ds.bound(i) = 0;
			if (sf.isValid())
				sf.set(ptoff, zero);
			if (smd.isValid())
				smd.set(ptoff, 0);
		}
		if (winklerp.isValid())
		{
			if (winklerp.get(ptoff) == 1)
				ds.winkler(i) = 1;
			else
				ds.winkler(i) = 0;
		}
		if (corep.isValid())
		{
			if (corep.get(ptoff) == 1)
				ds.core(i) = 1;
			else
				ds.core(i) = 0;

		}
	}
}

void gridUpdate::get_surface_force(const GU_Detail *gdp, const GA_Offset &primoff, const UT_Vector4i &tetra_points, Data_struct &ds)
{
	GA_ROHandleV3D force_handle(gdp->findAttribute(GA_ATTRIB_POINT, "surface_force"));
	GA_ROHandleI surface_handle(gdp->findAttribute(GA_ATTRIB_POINT, "surface"));
	unsigned int point_range = gdp->getPointRange().getEntries();
	if (force_handle.isValid() && surface_handle.isValid())
	{
		const GA_Primitive *prim = gdp->getPrimitive(primoff);
		int count = 0;
	}
}

void gridUpdate::global_indices(GU_Detail *gdp, Data_struct &ds)
{
	unsigned int point_count = gdp->getNumPoints();
	ds.global_index(0) = ds.tetra_points(0);
	ds.global_index(1) = ds.tetra_points(0) + point_count;
	ds.global_index(2) = ds.tetra_points(0) + 2 * point_count;
	ds.global_index(3) = ds.tetra_points(1);
	ds.global_index(4) = ds.tetra_points(1) + point_count;
	ds.global_index(5) = ds.tetra_points(1) + 2 * point_count;
	ds.global_index(6) = ds.tetra_points(2);
	ds.global_index(7) = ds.tetra_points(2) + point_count;
	ds.global_index(8) = ds.tetra_points(2) + 2 * point_count;
	ds.global_index(9) = ds.tetra_points(3);
	ds.global_index(10) = ds.tetra_points(3) + point_count;
	ds.global_index(11) = ds.tetra_points(3) + 2 * point_count;

	ds.global_bound_index(0) = ds.tetra_points(0);
	ds.global_bound_index(1) = ds.tetra_points(0);
	ds.global_bound_index(2) = ds.tetra_points(0);
	ds.global_bound_index(3) = ds.tetra_points(1);
	ds.global_bound_index(4) = ds.tetra_points(1);
	ds.global_bound_index(5) = ds.tetra_points(1);
	ds.global_bound_index(6) = ds.tetra_points(2);
	ds.global_bound_index(7) = ds.tetra_points(2);
	ds.global_bound_index(8) = ds.tetra_points(2);
	ds.global_bound_index(9) = ds.tetra_points(3);
	ds.global_bound_index(10) = ds.tetra_points(3);
	ds.global_bound_index(11) = ds.tetra_points(3);
}


void gridUpdate::assemble_global_implicit(GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds, SIM_Data_Sparse *sparse, const Material &parms, const Wilson &wil, const bool &use_viscosity)
{
	GA_RWHandleM3D stress_handle(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "cauchy_stress"));
	GA_ROHandleV3D force_handle(gdp->findAttribute(GA_ATTRIB_POINT, "volume_force"));
	GA_ROHandleV3D lp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P0"));
	GA_ROHandleV3D pp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P1"));
	GA_ROHandleV3D p_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
	GA_ROHandleD smd(gdp->findAttribute(GA_ATTRIB_POINT, "surface_mass_density"));
	GA_ROHandleV3D vp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "v"));
	GA_ROHandleV3D sf(gdp->findAttribute(GA_ATTRIB_POINT, "surface_force"));
	GA_ROHandleD dc(gdp->findAttribute(GA_ATTRIB_POINT, "density_contrast"));
	GA_ROHandleV3D displacement_p(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));
	GA_ROHandleD gp_element(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "g0"));
	GA_ROHandleV3 disp(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
	GA_ROHandleD element_density(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "density"));
	GA_RWHandleV3 winkler_test(gdp->findAttribute(GA_ATTRIB_POINT, "wt"));
	GA_RWHandleV3D gvp(gdp->findAttribute(GA_ATTRIB_POINT, "gravity_vector"));
	GA_RWHandleI core_p(gdp->findAttribute(GA_ATTRIB_POINT, "core"));

	if (stress_handle.isValid() && vp_handle.isValid() && pp_handle.isValid() && lp_handle.isValid())
	{
		unsigned int point_range = gdp->getPointRange().getEntries();

		double prim_surf_area = 0;
		double a = ds.determinant *1.0/ 120.0;
		double b = ds.determinant * 1 / 6.0;
		double c = ds.surf_area / 12.0;

		UT_Vector3D prim_normal(0.0);
		if (ds.winkler(0) + ds.winkler(1) + ds.winkler(2)  == 3)
		{
			UT_Vector3D d0 = gdp->getPos3(ds.tetra_points_offset(0)) - gdp->getPos3(ds.tetra_points_offset(1));
			UT_Vector3D d1 = gdp->getPos3(ds.tetra_points_offset(0)) - gdp->getPos3(ds.tetra_points_offset(2));

			prim_normal = cross(d0, d1);
			prim_normal.normalize();
		}
		if (ds.winkler(3) + ds.winkler(1) + ds.winkler(2) == 3)
		{
			UT_Vector3D d0 = gdp->getPos3(ds.tetra_points_offset(3)) - gdp->getPos3(ds.tetra_points_offset(1));
			UT_Vector3D d1 = gdp->getPos3(ds.tetra_points_offset(3)) - gdp->getPos3(ds.tetra_points_offset(2));

			prim_normal = cross(d0, d1);
			prim_normal.normalize();
		}
		if (ds.winkler(0) + ds.winkler(3) + ds.winkler(2) == 3)
		{
			UT_Vector3D d0 = gdp->getPos3(ds.tetra_points_offset(0)) - gdp->getPos3(ds.tetra_points_offset(3));
			UT_Vector3D d1 = gdp->getPos3(ds.tetra_points_offset(0)) - gdp->getPos3(ds.tetra_points_offset(2));

			prim_normal = cross(d0, d1);
			prim_normal.normalize();
		}
		if (ds.winkler(0) + ds.winkler(1) + ds.winkler(3) == 3)
		{
			UT_Vector3D d0 = gdp->getPos3(ds.tetra_points_offset(0)) - gdp->getPos3(ds.tetra_points_offset(1));
			UT_Vector3D d1 = gdp->getPos3(ds.tetra_points_offset(0)) - gdp->getPos3(ds.tetra_points_offset(3));

			prim_normal = cross(d0, d1);
			prim_normal.normalize();
		}

		GA_Primitive *prim = gdp->getPrimitive(primoff);
	
		unsigned int i = 0;
		unsigned int bi = 0;
		
		UT_Vector strain_vector(0, 5);
		ds.stress = stress_handle.get(primoff);
		ds.B.postMult(ds.point_values, strain_vector);
		ds.M2.postMult(strain_vector, ds.stress_vector);
	
		add_full_stress_vector(ds, gdp, prim, use_viscosity);
		form_stress_tensor(ds);
		stress_handle.set(primoff, ds.stress);
		ds.point_values.zero();
		ds.BT.postMult(ds.stress_vector, ds.point_values);

		ds.BT.postMult(ds.N, ds.pre_stress);
		UT_MatrixD pre_stress_t(0, 11, 0, 11);

		ds.pre_stress.transpose(pre_stress_t);
		ds.pre_stress.addScaledMatrix(pre_stress_t, 1);

		if (!ds.stashed)
		{
			UT_Matrix diagonal(0, 11, 0, 11);
			diagonal.zero();
			assemble_linear_element_matrix(ds, use_viscosity);	
			ds.non_linear_element_3D = ds.element_stiffness_matrix_nl;
			ds.non_linear_element_3D.addScaledMatrix(ds.element_stiffness_matrix, 1);
			
			for (unsigned int i = 0; i < 4; i++)
			{
				rotate_degree_of_freedom(gdp, primoff, ds, i, ds.stashed);		
				add_winkler(gdp, c, prim_normal, ds, i);
			}
			
		}
		
		UT_Vector vector_prestress(0, 11);
		vector_prestress.zero();

		UT_Vector vd(0, 11);
		vd.zero();

		ds.pre_stress.postMult(ds.vector_displacements, vector_prestress);
		vector_prestress *= abs(element_density(primoff));
		ds.pre_stress.postMult(ds.divergence_vector, vd);

		unsigned int q = 0;
		unsigned int i_1 = 0;
		const double mass = 5 * parms.density;

		for (GA_Iterator pointit(prim->getPointRange()); !pointit.atEnd(); ++pointit)
		{
			UT_Vector3D force_vector(0.0);
			UT_Vector3D surface_force_vector(0.0);
			double surface_mass_density(0.0);
			
			UT_Vector3D posi = p_handle(ds.tetra_points(i_1));
			double Li = posi.length();
			posi = posi / Li;

			if (ds.core(0) + ds.core(1) + ds.core(2) + ds.core(3) == 3)
			{
				calculate_and_apply_winkler(gdp, surface_force_vector, ds, i_1);
				UT_Vector3D n_p = posi / posi.length();
				double scaled_c = c*abs(dot(n_p, prim_normal));
				surface_force_vector *= 1 / c*scaled_c;
			}

			ds.lumped_mass(ds.tetra_points(i_1)) += a*mass;
	
			UT_Vector3D pos0 = p_handle(ds.tetra_points(0));
			double L0 = pos0.length();
			pos0 = pos0 / L0;

			UT_Vector3D pos1 = p_handle(ds.tetra_points(1));
			double L1 = pos1.length();
			pos1 = pos1 / L1;

			UT_Vector3D pos2 = p_handle(ds.tetra_points(2));
			double L2 = pos2.length();
			pos2 = pos2 / L2;

			UT_Vector3D pos3 = p_handle(ds.tetra_points(3));
			double L3 = pos3.length();
			pos3 = pos3 / L3;

			UT_Vector3D Bui = pos0*(vd(0) + vd(0 + 1) + vd(0 + 2));
			UT_Vector3D Bu0 = pos1*(vd(3) + vd(3 + 1) + vd(3 + 2));
			UT_Vector3D Bu1 = pos2*(vd(6) + vd(6 + 1) + vd(6 + 2));
			UT_Vector3D Bu2 = pos3*(vd(9) + vd(9 + 1) + vd(9 + 2));
			UT_Vector3D Bu3 = posi*(vd(3*i_1) + vd(3*i_1 + 1) + vd(3*i_1 + 2));

			if (sf.isValid() && ds.bound(0) + ds.bound(1) + ds.bound(2) + ds.bound(3) == 3)
				surface_force_vector += ds.bound(i_1)*(sf(ds.tetra_points_offset(i_1)) + sf(ds.tetra_points_offset(0)) + sf(ds.tetra_points_offset(1)) + sf(ds.tetra_points_offset(2)) + sf(ds.tetra_points_offset(3)));
					
			if (force_handle.isValid())
				force_vector += force_handle(ds.tetra_points_offset(i_1)) + force_handle(ds.tetra_points_offset(0)) + force_handle(ds.tetra_points_offset(1)) + force_handle(ds.tetra_points_offset(2)) + force_handle(ds.tetra_points_offset(3));
			
			UT_Vector3 external_force(0.0);
			UT_Vector3 winkler_attp(0.0);
		
			UT_Vector3D surface_vector = p_handle(ds.tetra_points_offset(i_1));
			UT_Vector3D surface_normal = surface_vector / surface_vector.length();
			
			double L = p_handle(ds.tetra_points_offset(i_1)).length();

			for (unsigned int j = 0; j < 3; j++)
			{
				if (disp.isInvalid())
				{
					ds.global_force(ds.global_index(q)) -= b*(ds.point_values(3 * i_1 + j));
					ds.global_force(ds.global_index(q)) -= c*(surface_force_vector(j));
					ds.global_force(ds.global_index(q)) -= abs(b)*(vector_prestress(3 * i_1 + j));
					ds.global_force(ds.global_index(q)) -= gp_element(primoff)*element_density(primoff)*abs(a)*(Bui(j) + Bu0(j) + Bu1(j) + Bu2(j) + Bu3(j));
					//winkler_attp(j) = -abs(b)*vector_prestress(3 * i_1 + j);
				}
				else
				{ 
					const UT_Vector3 boundary_conditions = disp(ds.tetra_points_offset(i_1));
					if (boundary_conditions(j) == 0)
					{
						ds.global_force(ds.global_index(q)) -= b*(ds.point_values(3 * i_1 + j));
						ds.global_force(ds.global_index(q)) -= c*(surface_force_vector(j));
						ds.global_force(ds.global_index(q)) -= abs(b)*(vector_prestress(3 * i_1 + j));
						ds.global_force(ds.global_index(q)) -= gp_element(primoff)*element_density(primoff)*abs(a)*(Bui(j) + Bu0(j) + Bu1(j) + Bu2(j) + Bu3(j));
					//	winkler_attp(j) = -abs(b)*vector_prestress(3 * i_1 + j);
					}
				}

				if (!ds.stashed)
				{ 
					for (unsigned int w = 0; w < 12; w++)
					{
						inject_element_matrix(gdp, w, q, ds, sparse, wil);
					}
				}
				q += 1;
			}
		
			i_1 += 1;
		}
	}
}

void gridUpdate::assemble_global( GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds, SIM_Data_Sparse *sparse, const Material &parms, const bool &use_viscosity)
{
	GA_RWHandleM3D stress_handle(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "cauchy_stress"));
	GA_ROHandleV3D force_handle(gdp->findAttribute(GA_ATTRIB_POINT, "volume_force"));
	GA_ROHandleV3D sf(gdp->findAttribute(GA_ATTRIB_POINT, "surface_force"));
	GA_ROHandleV3D lp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P0"));
	GA_ROHandleV3D pp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P1"));
	GA_ROHandleV3D vp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "v"));
	
	if (stress_handle.isValid() && vp_handle.isValid() && pp_handle.isValid() && lp_handle.isValid())
	{
		unsigned int point_range = gdp->getPointRange().getEntries();

		double prim_surf_area = 0;
		double a = ds.determinant / 120.0;
		double b = ds.determinant / 6.0;
		double c = ds.surf_area / 12.0;


		GA_Primitive *prim = gdp->getPrimitive(primoff);
		unsigned int i = 0;
		UT_Vector3i boundary_points;
		unsigned int boundary_check = 0;
		UT_Vector point_values(0, 11);


		for (GA_Iterator pointit(prim->getPointRange()); !pointit.atEnd(); ++pointit)
		{
			GA_Offset ptoff = *pointit;
			bool is_b = Collider::is_surface(gdp, ptoff, "surface_points");
			if (is_b && i < 3)
			{
				boundary_check += 1;
				boundary_points(i) = ptoff;
			}

			for (unsigned int j = 0; j < 3; j++)
			{
				point_values(3 * i + j) = pp_handle(ds.tetra_points_offset(i))(j) - lp_handle(ds.tetra_points_offset(i))(j);
			}
			i += 1;
		}
		UT_Vector strain_vector(0, 5);
		ds.stress = stress_handle.get(primoff);
		ds.B.postMult(point_values, strain_vector);
		ds.sigma.postMult(strain_vector, ds.stress_vector);
		add_full_stress_vector(ds, gdp, prim, use_viscosity);
		form_stress_tensor(ds);
		stress_handle.set(primoff, ds.stress);

		point_values.zero();
		ds.BT.postMult(ds.stress_vector, point_values);

		unsigned int i_1 = 0;
		for (GA_Iterator pointit(prim->getPointRange()); !pointit.atEnd(); ++pointit)
		{
			GA_Offset ptoff = *pointit;
			GA_Index ptindex = pointit.getIndex();

			double mass = (5 * parms.density);
			ds.lumped_mass(ptindex) += a*mass;

			UT_Vector3D force_vector(0.0);
			UT_Vector3D surface_force_vector(0.0);

			if (sf.isValid())
			{
				surface_force_vector = ds.bound(i_1)*(sf(ds.tetra_points_offset(i_1)) + sf(ds.tetra_points_offset(0)) + sf(ds.tetra_points_offset(1)) + sf(ds.tetra_points_offset(2)) + sf(ds.tetra_points_offset(3)));
			}
			if (force_handle.isValid())
			{
				force_vector = force_handle(ds.tetra_points_offset(i_1)) + force_handle(ds.tetra_points_offset(0)) + force_handle(ds.tetra_points_offset(1)) + force_handle(ds.tetra_points_offset(2)) + force_handle(ds.tetra_points_offset(3));
			}

			UT_Vector3 external_force(0.0);
			
			for (unsigned int j = 0; j < 3; j++)
			{
				ds.global_force(ptindex + j*(point_range)) += a*(force_vector(j) + external_force(j));
				ds.global_force(ptindex + j*(point_range)) -= b*point_values(3 * i_1 + j);
				ds.global_force(ptindex + j*(point_range)) -= c*(surface_force_vector(j));

			}
			i_1 += 1;
		}
	}
}
void gridUpdate::update_jacobians(GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds)
{
	GA_ROHandleV3D pp(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
	GA_RWHandleD detp(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "determinant"));
	GA_RWHandleM4D invjp(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "inv_j"));
	GA_RWHandleD surfp(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "surface_area"));

	if (detp.isValid() && invjp.isValid() && surfp.isValid())
	{

		ds.determinant = detp(primoff);
		ds.surf_area = surfp(primoff);
		ds.inverse_jacobian = invjp(primoff);

	}
	if(!abs(ds.determinant) || !ds.stashed)
	{
		UT_Vector3D p0 = pp.get(ds.tetra_points_offset(0));
		UT_Vector3D p1 = pp.get(ds.tetra_points_offset(1));
		UT_Vector3D p2 = pp.get(ds.tetra_points_offset(2));
		UT_Vector3D p3 = pp.get(ds.tetra_points_offset(3));

		ds.jacobian(0, 0) = 1;
		ds.jacobian(1, 0) = p0(0);
		ds.jacobian(2, 0) = p0(1);
		ds.jacobian(3, 0) = p0(2);

		ds.jacobian(0, 1) = 1;
		ds.jacobian(1, 1) = p1(0);
		ds.jacobian(2, 1) = p1(1);
		ds.jacobian(3, 1) = p1(2);

		ds.jacobian(0, 2) = 1;
		ds.jacobian(1, 2) = p2(0);
		ds.jacobian(2, 2) = p2(1);
		ds.jacobian(3, 2) = p2(2);

		ds.jacobian(0, 3) = 1;
		ds.jacobian(1, 3) = p3(0);
		ds.jacobian(2, 3) = p3(1);
		ds.jacobian(3, 3) = p3(2);

		ds.determinant = -abs(ds.jacobian.determinant());
		ds.inverse_jacobian = ds.jacobian;
		ds.inverse_jacobian.invertDouble();


		if (!ds.bound(0) && ds.bound(1) && ds.bound(2) && ds.bound(3))
			ds.surf_area = 0.5*cross(p3 - p1, p2 - p1).length();
		else if (!ds.bound(1) && ds.bound(0) && ds.bound(2) && ds.bound(3))
			ds.surf_area = 0.5*cross(p3 - p0, p2 - p0).length();
		else if (!ds.bound(2) && ds.bound(1) && ds.bound(0) && ds.bound(3))
			ds.surf_area = 0.5*cross(p3 - p0, p1 - p0).length();
		else if (!ds.bound(3) && ds.bound(1) && ds.bound(2) && ds.bound(0))
			ds.surf_area = 0.5*cross(p1 - p0, p2 - p0).length();
	

		if (!ds.winkler(0) && ds.winkler(1) && ds.winkler(2) && ds.winkler(3))
		{
			UT_Vector3D normal = 0.5*cross(p3 - p1, p2 - p1);
			ds.surf_area = normal.length();
		}
		else if (!ds.winkler(1) && ds.winkler(0) && ds.winkler(2) && ds.winkler(3))
		{
			UT_Vector3D normal = 0.5*cross(p3 - p0, p2 - p0);
			ds.surf_area = normal.length();
		}
		else if (!ds.winkler(2) && ds.winkler(1) && ds.winkler(0) && ds.winkler(3))
		{
			UT_Vector3D normal = 0.5*cross(p3 - p0, p1 - p0);
			ds.surf_area = normal.length();
		}
		else if (!ds.winkler(3) && ds.winkler(1) && ds.winkler(2) && ds.winkler(0))
		{
			UT_Vector3D normal = 0.5*cross(p1 - p0, p2 - p0);
			ds.surf_area = normal.length();	
		}
		if (detp.isValid())
			detp.set(primoff, ds.determinant);
		if (surfp.isValid())
			surfp.set(primoff, ds.surf_area);
		if (invjp.isValid())
			invjp.set(primoff, ds.inverse_jacobian);
	}

}

void gridUpdate::assemble_element_matrices(GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds, const Material &parms, int __explicit__)
{
	GA_ROHandleV3D lp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P0"));
	GA_ROHandleV3D pp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P1"));
	GA_ROHandleV3D p0_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
	GA_ROHandleV3D d_handle(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));
	GA_RWHandleV3D wi(gdp->findAttribute(GA_ATTRIB_POINT, "wt"));
	GA_ROHandleD gp_element(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "g0"));
	GA_ROHandleD element_density(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "density"));
	GA_RWHandleD potential_p(gdp->findAttribute(GA_ATTRIB_POINT, "potential"));

	ds.M1.zero();

	unsigned int loop_max = 6;

	unsigned int bound = loop_max / 2;
	double tau =  2*parms.viscosity*(1.0 + parms.poisson_ratio) / parms.youngs_modulus;

	const double aa = (1 - exp(-parms.dt / tau))  * (2.0 / 3.0);
	const double ab = (1 - exp(-parms.dt / tau)) * (-1.0 / 3.0);
	const double ac = (1 - exp(-parms.dt / tau));

	ds.M2.zero();

	const double ah = parms.youngs_modulus / (3.0 * (1.0 - 2.0 * parms.poisson_ratio));
	const double bh = parms.youngs_modulus / (3.0 * (1.0 - 2.0 * parms.poisson_ratio));

	const double af = 4.0 * parms.youngs_modulus / (6.0 * (1.0 + parms.poisson_ratio)) * (tau / parms.dt)*(1.0 - exp(-parms.dt / tau));
	const double bf = -2.0 * parms.youngs_modulus / (6.0 * (1.0 + parms.poisson_ratio)) * (tau / parms.dt)*(1.0 - exp(-parms.dt / tau));
	const double cf = 3.0 * parms.youngs_modulus / (6.0 * (1.0 + parms.poisson_ratio)) * (tau / parms.dt)*(1.0 - exp(-parms.dt / tau));

	ds.sigma.zero();

	const double diagonal1 = 1 - parms.poisson_ratio;
	const double diagonal2 = 0.5 - parms.poisson_ratio;
	const double diagonal3 = parms.poisson_ratio;
	const double multiplier = parms.youngs_modulus / ((1 + parms.poisson_ratio)*(1 - 2 * parms.poisson_ratio));

	const double diag1 = diagonal1 * multiplier;
	const double diag2 = diagonal2 * multiplier;
	const double diag3 = diagonal3 * multiplier;
	GA_Offset primoff = prim->getMapOffset();
	if (__explicit__ == -2)
		loop_max = 4;

	for (unsigned int i = 0; i < loop_max; i++)
	{
		gridUpdate::update_tetrapoints(gdp, prim, ds, i);
	
			for (unsigned int j = 0; j < loop_max; j++)
			{
				if (i < 4 && j < 3)
				{
					ds.point_values(3 * i + j) = pp_handle(ds.tetra_points_offset(i))(j) - lp_handle(ds.tetra_points_offset(i))(j);
					if (gp_element.isValid() && element_density.isValid() && d_handle.isValid() && pp_handle.isValid())
					{
						UT_Vector3D pos = pp_handle(ds.tetra_points_offset(i));
						UT_Vector3D di = d_handle(ds.tetra_points_offset(i));
						fpreal L = pos.length();
						UT_Vector3D pos_n = -pos / L;
							
						fpreal sn = pos_n.dot(di);
						ds.vector_displacements(3 * i + j) = abs(gp_element(primoff))*sn;
						if(potential_p.isValid())
							ds.vector_displacements(3 * i + j) += potential_p(ds.tetra_points_offset(i));
						
						ds.divergence_vector(3 * i + j) = di(j);
					}
				}
			if (__explicit__ != -2)
			{
				if (i == j && i < bound)
				{
					ds.M1(i, j) = 1.0 - aa;
					ds.M2(i, j) = af + ah;
					ds.sigma(i, j) = diag1;
				}
				if (i == j && i >= bound)
				{
					ds.M1(i, j) = 1.0 - ac;
					ds.M2(i, j) = cf;
					ds.sigma(i, j) = diag2;
				}
				if (i != j && i < bound && j < bound)
				{
					ds.M1(i, j) = -ab;
					ds.M2(i, j) = bf + bh;
					ds.sigma(i, j) = diag3;
				}
			}
		}
	}
}

void gridUpdate::assemble_element_matrices_laplacian(GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds, const Material &parms)
{
	GA_ROHandleV3D d_handle(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));

	unsigned int loop_max = 4;
	unsigned int bound = loop_max / 2;
	GA_Offset primoff = prim->getMapOffset();

	for (unsigned int i = 0; i < loop_max; i++)
	{
		gridUpdate::update_tetrapoints(gdp, prim, ds, i);

		for (unsigned int j = 0; j < loop_max; j++)
		{
			if (i < 4 && j < 3)
			{
				if (d_handle.isValid())
				{

					UT_Vector3D di = d_handle(ds.tetra_points_offset(i));
					ds.divergence_vector(3 * i + j) = di(j);
				}
			}
		}
	}
}


void gridUpdate::assemble_m1(Data_struct &ds, const Material &parms)
{
	ds.M1.zero();
	unsigned int size = ds.sigma.columns();
	unsigned int loop_max = 6;

	unsigned int bound = loop_max / 2;
	float tau = 2 * parms.viscosity*(1 + parms.poisson_ratio) / parms.youngs_modulus;

	float aa = (1 - exp(-parms.dt / tau))  * (2.0 / 3.0);
	float ab = (1 - exp(-parms.dt / tau)) * (-1.0 / 3.0);
	float ac = (1 - exp(-parms.dt / tau)) * (1.0);

	for (unsigned int i = 0; i < loop_max; i++)
	{
		for (unsigned int j = 0; j < loop_max; j++)
		{
			if (i == j && i < bound)
				ds.M1(i, j) = 1.0 - aa;
			if (i == j && i >= bound)
				ds.M1(i, j) = 1.0 - ac;
			if (i != j && i < bound && j < bound)
				ds.M1(i, j) = -ab;
		}
	}
};
void gridUpdate::assemble_m2(Data_struct &ds, const Material &parms)
{
	ds.M2.zero();
	unsigned int size = ds.sigma.columns();
	unsigned int loop_max = 6;
	
	unsigned int bound = loop_max / 2;
	float ah = parms.youngs_modulus / (3.0 * (1.0 - 2.0 * parms.poisson_ratio));
	float bh = parms.youngs_modulus / (3.0 * (1.0 - 2.0 * parms.poisson_ratio));

	double tau = 2.0 * parms.viscosity*(1 + parms.poisson_ratio) / parms.youngs_modulus;

	float af = 4.0 * parms.youngs_modulus / (6.0 * (1.0 + parms.poisson_ratio)) * (tau / parms.dt)*(1.0 - exp(-parms.dt / tau));
	float bf = -2.0 * parms.youngs_modulus / (6.0 * (1.0 + parms.poisson_ratio)) * (tau / parms.dt)*(1.0 - exp(-parms.dt / tau));
	float cf = 3.0 * parms.youngs_modulus / (6.0 * (1.0 + parms.poisson_ratio)) * (tau / parms.dt)*(1.0 - exp(-parms.dt / tau));

	for (unsigned int i = 0; i < loop_max; i++)
	{
		for (unsigned int j = 0; j < loop_max; j++)
		{
			if (i == j && i < bound)
				ds.M2(i, j) = af + ah;
			if (i == j && i >= bound)
				ds.M2(i, j) = cf;
			if (i != j && i < bound && j < bound)
				ds.M2(i, j) = bf + bh;
		}
	}
};
void gridUpdate::assemble_E(Data_struct &ds, const Material &parms)
{
	ds.sigma.zero();
	unsigned int size = ds.sigma.columns();
	unsigned int loop_max = 6;
	if (size < 6)
		loop_max = 4;

	unsigned int bound = loop_max / 2;

	float diagonal1 = 1 - parms.poisson_ratio;
	float diagonal2 = 0.5 - parms.poisson_ratio;
	float diagonal3 = parms.poisson_ratio;
	float multiplier = parms.youngs_modulus / ((1 + parms.poisson_ratio)*(1 - 2* parms.poisson_ratio));

	float diag1 = diagonal1 * multiplier;
	float diag2 = diagonal2 * multiplier;
	float diag3 = diagonal3 * multiplier;

	for (unsigned int i = 0; i < loop_max; i++)
	{
		for (unsigned int j = 0; j < loop_max; j++)
		{
			if (i == j && i < bound)
				ds.sigma(i, j) = diag1;
			if (i == j && i >= bound)
				ds.sigma(i, j) = diag2;
			if (i != j && i < bound && j < bound)
				ds.sigma(i, j) = diag3;
		}
	}
};
void gridUpdate::update_material_parms(const GU_Detail *gdp, const GA_Primitive *prim, Data_struct &ds, Material &parms, bool &attribute_present, bool &use_vic)
{

	GA_ROHandleD densityp(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "density"));
	GA_ROHandleD viscosityp(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "viscosity"));
	GA_ROHandleD poissonp(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "poisson"));
	GA_ROHandleD youngp(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "Y"));
	GA_Offset primoff = prim->getMapOffset();

	if (parms.viscosity < 0.0000001)
		use_vic = false;
	else
		use_vic = true;

	if (viscosityp.isValid())
		parms.viscosity = abs(viscosityp(primoff)), attribute_present = true;

	if (poissonp.isValid())
		parms.poisson_ratio = abs(poissonp(primoff)), attribute_present = true;

	if (youngp.isValid())
		parms.youngs_modulus = abs(youngp(primoff)), attribute_present = true;

	if (densityp.isValid())
		parms.density = abs(densityp(primoff)), attribute_present = true;

	if (parms.poisson_ratio <= 0)
		parms.poisson_ratio = 0.001;
}
void gridUpdate::add_full_stress_vector(Data_struct &ds, GU_Detail *gdp, const GA_Primitive *prim, const bool &use_vic)
{
	UT_Vector lingering_stress(0,5);
	lingering_stress.zero();
	lingering_stress(0) = ds.stress(0, 0);
	lingering_stress(1) = ds.stress(1, 1);
	lingering_stress(2) = ds.stress(2, 2);

	lingering_stress(3) = ds.stress(1, 2);
	lingering_stress(4) = ds.stress(2, 0);
	lingering_stress(5) = ds.stress(0, 1);

	if (use_vic)
	{
		ds.M1.postMult(lingering_stress, ds.m1_holder);
		ds.stress_vector += ds.m1_holder;
	}
	
}

void gridUpdate::form_stress_tensor(Data_struct &ds)
{
	if (ds.stress_vector.length() == 6)
	{
		ds.stress(0, 0) = ds.stress_vector(0);
		ds.stress(1, 1) = ds.stress_vector(1);
		ds.stress(2, 2) = ds.stress_vector(2);

		ds.stress(0, 2) = ds.stress_vector(4);
		ds.stress(2, 0) = ds.stress_vector(4);

		ds.stress(1, 0) = ds.stress_vector(5);
		ds.stress(0, 1) = ds.stress_vector(5);

		ds.stress(2, 1) = ds.stress_vector(3);
		ds.stress(1, 2) = ds.stress_vector(3);
	}
}
void gridUpdate::assemble_linear_element_matrix(Data_struct &ds, bool use_vic)
{
	UT_MatrixD B_temp;
	B_temp.init(0, 5, 0, 11);
	if (!use_vic)
	{
		ds.sigma.postMult(ds.B, B_temp);
	}
	else
	{
		ds.M2.postMult(ds.B, B_temp);
	}
	ds.element_stiffness_matrix.zero();
	ds.BT.postMult(B_temp, ds.element_stiffness_matrix);
}

void gridUpdate::inject_element_matrix(const GU_Detail *gdp, unsigned int i, unsigned int j, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &wil)
{
	GA_ROHandleV3 pin_p(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
	GA_ROHandleD dc(gdp->findAttribute(GA_ATTRIB_POINT, "density_contrast"));
	GA_ROHandleI cc(gdp->findAttribute(GA_ATTRIB_POINT, "core"));
	GA_ROHandleV3D pp(gdp->findAttribute(GA_ATTRIB_POINT, "P"));

	unsigned int point_count = gdp->getPointRange().getEntries();
	fpreal value = 0;
	fpreal pre_stress_value = 0;
	
	unsigned int x = ds.global_index(i);
	unsigned int y = ds.global_index(j);

	unsigned int offset1 = ds.global_bound_index(i);
	unsigned int offset2 = ds.global_bound_index(j);
	
	const GA_PointGroup *found_group = gdp->findPointGroup("surface_points");
	value = (ds.non_linear_element_3D(i, j));
	value *= ds.determinant * 1.0 / 6.0;
	
	sparse->stiffness_matrix->addToElement(x, y, value);
}
double gridUpdate::point_to_prim(const GA_Primitive *prim, GA_ROHandleD &attrib_handle, const Data_struct &ds)
{
	fpreal prim_attrib = 0;

	for (GA_Iterator pointit(prim->getPointRange()); !pointit.atEnd(); ++pointit)
	{
		GA_Offset ptoff = *pointit;
		prim_attrib += attrib_handle.get(ptoff);
	}
	prim_attrib *= 0.25;
	return prim_attrib;
}

void gridUpdate::assemble_B(Data_struct &ds)
{
	ds.B.zero();

	ds.B(0, 0) = ds.inverse_jacobian(0, 1);
	ds.B(0, 3) = ds.inverse_jacobian(1, 1);
	ds.B(0, 6) = ds.inverse_jacobian(2, 1);
	ds.B(0, 9) = ds.inverse_jacobian(3, 1);

	ds.B(1, 1 + 0) = ds.inverse_jacobian(0, 2);
	ds.B(1, 1 + 3) = ds.inverse_jacobian(1, 2);
	ds.B(1, 1 + 6) = ds.inverse_jacobian(2, 2);
	ds.B(1, 1 + 9) = ds.inverse_jacobian(3, 2);

	ds.B(2, 2 + 0) = ds.inverse_jacobian(0, 3);
	ds.B(2, 2 + 3) = ds.inverse_jacobian(1, 3);
	ds.B(2, 2 + 6) = ds.inverse_jacobian(2, 3);
	ds.B(2, 2 + 9) = ds.inverse_jacobian(3, 3);

	ds.B(3, 1 + 0) = ds.inverse_jacobian(0, 1);
	ds.B(3, 1 + 3) = ds.inverse_jacobian(1, 1);
	ds.B(3, 1 + 6) = ds.inverse_jacobian(2, 1);
	ds.B(3, 1 + 9) = ds.inverse_jacobian(3, 1);

	ds.B(3, 0) = ds.inverse_jacobian(0, 2);
	ds.B(3, 3) = ds.inverse_jacobian(1, 2);
	ds.B(3, 6) = ds.inverse_jacobian(2, 2);
	ds.B(3, 9) = ds.inverse_jacobian(3, 2);

	ds.B(4, 1 + 0) = ds.inverse_jacobian(0, 3);
	ds.B(4, 1 + 3) = ds.inverse_jacobian(1, 3);
	ds.B(4, 1 + 6) = ds.inverse_jacobian(2, 3);
	ds.B(4, 1 + 9) = ds.inverse_jacobian(3, 3);

	ds.B(4, 2 + 0) = ds.inverse_jacobian(0, 2);
	ds.B(4, 2 + 3) = ds.inverse_jacobian(1, 2);
	ds.B(4, 2 + 6) = ds.inverse_jacobian(2, 2);
	ds.B(4, 2 + 9) = ds.inverse_jacobian(3, 2);

	ds.B(5, 2 + 0) = ds.inverse_jacobian(0, 1);
	ds.B(5, 2 + 3) = ds.inverse_jacobian(1, 1);
	ds.B(5, 2 + 6) = ds.inverse_jacobian(2, 1);
	ds.B(5, 2 + 9) = ds.inverse_jacobian(3, 1);

	ds.B(5, 0) = ds.inverse_jacobian(0, 3);
	ds.B(5, 3) = ds.inverse_jacobian(1, 3);
	ds.B(5, 6) = ds.inverse_jacobian(2, 3);
	ds.B(5, 9) = ds.inverse_jacobian(3, 3);

	ds.B.transpose(ds.BT);
	
};

void gridUpdate::apply_inertia(Data_struct &ds, SIM_Data_Sparse *sparse, GU_Detail *gdp, const Wilson &wil)
{
	GA_ROHandleV3 disp(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
	GA_ROHandleV3D vp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "v"));
	GA_ROHandleV3D accp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "acc"));
	GA_ROHandleI core(gdp->findAttribute(GA_ATTRIB_POINT, "core"));
	unsigned int pointCount = gdp->getNumPoints();
	
	
	for (GA_Iterator point_iter(gdp->getPointRange()); !point_iter.atEnd(); ++point_iter)
	{
		const GA_Offset offset = *point_iter;
		const GA_Index index = point_iter.getIndex();
		UT_Vector3D acc = accp_handle(offset);
		UT_Vector3D v = vp_handle(offset);
		bool	is_core = false;
		if (core.isValid())
		{
			if(core(offset) )
				is_core = true;
		}
		fpreal mass_contrib = 0;

		mass_contrib = wil.a0*ds.lumped_mass(index);
		UT_Vector3D inertia;

		inertia = ds.lumped_mass(index)*(wil.a1*v + wil.a2*acc);
		rotate_degree_of_freedom(gdp, ds, offset, index, ds.global_force, true);
		if (disp.isInvalid())
		{
			if (!ds.stashed)
			{
				sparse->stiffness_matrix->addToElement(index, index, mass_contrib);
				sparse->stiffness_matrix->addToElement(index + pointCount, index + pointCount, mass_contrib);
				sparse->stiffness_matrix->addToElement(index + 2 * pointCount, index + 2 * pointCount, mass_contrib);
			}
			ds.fixed_global_force(index) = inertia(0);
			ds.fixed_global_force(index + pointCount) = inertia(1);
			ds.fixed_global_force(index + 2 * pointCount) = inertia(2);
		}
		
		if (disp.isValid())
		{
			UT_Vector3 di = disp(offset);
			if (di(0) != 1)
			{
				if(!ds.stashed)
					sparse->stiffness_matrix->addToElement(index, index, mass_contrib);
			}
			if (di(1) != 1)
			{
				if(!ds.stashed)
					sparse->stiffness_matrix->addToElement(index + pointCount, index + pointCount, mass_contrib);
			}
			if (di(2) != 1)
			{
				if(!ds.stashed)
					sparse->stiffness_matrix->addToElement(index + 2 * pointCount, index + 2 * pointCount, mass_contrib);
			}

			if (di(0) == 1 && !is_core)
			{
				if(!ds.stashed)
					sparse->stiffness_matrix->addToElement(index, index, 1);
				ds.global_force(index) = 0;
			}

			if (di(1) == 1 && !is_core)
			{
				if(!ds.stashed)
					sparse->stiffness_matrix->addToElement(index + pointCount, index + pointCount, 1);
				ds.global_force(index + 1 * pointCount) = 0;
			}

			if (di(2) == 1 && !is_core)
			{
				if(!ds.stashed)
					sparse->stiffness_matrix->addToElement(index + 2 * pointCount, index + 2 * pointCount, 1);
				ds.global_force(index + 2 * pointCount) = 0;
			}
		}
	}
}

void gridUpdate::step_solution_exp(GU_Detail *gdp, Data_struct &ds, Material &parms, const Wilson &wil)
{
	fpreal t_1 = pow(parms.dt, 2.0);
	unsigned int points = gdp->getNumPoints();
	GA_RWHandleV3D pp0(gdp->findAttribute(GA_ATTRIB_POINT, "P0"));
	GA_RWHandleV3D pp(gdp->findAttribute(GA_ATTRIB_POINT, "P1"));
	GA_RWHandleV3D vp(gdp->findAttribute(GA_ATTRIB_POINT, "v"));
	GA_RWHandleV3D accp(gdp->findAttribute(GA_ATTRIB_POINT, "acc"));
	GA_RWHandleV3D disp(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));
	GA_RWHandleI pin(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));

	if (pp0.isValid() && pp.isValid() && accp.isValid() && vp.isValid())
	{
		for (GA_Iterator pointit(gdp->getPointRange()); !pointit.atEnd(); ++pointit)
		{
			GA_Offset ptoff = *pointit;
			GA_Index index = pointit.getIndex();

			UT_Vector3D point_force(ds.global_force(index), ds.global_force(index + points), ds.global_force(index + 2 * points));
			UT_Vector3D p_next = 1.0 / ds.lumped_mass(index)*t_1*point_force + 2 * pp(ptoff) - pp0(ptoff);

			UT_Vector3D new_acc = (wil.a3*(p_next - pp(ptoff)) + wil.a4*vp(ptoff) + wil.a5*accp(ptoff));
			UT_Vector3D new_v = vp(ptoff) + wil.a6*(new_acc + accp(ptoff));

			vp.set(ptoff, new_v);
			accp.set(ptoff, new_acc);
			pp0.set(ptoff, pp(ptoff));
			if (disp.isValid())
				disp.add(ptoff, p_next - pp(ptoff));

			if (pin.isValid())
			{
				if (!pin(ptoff))
					pp.set(ptoff, p_next);
			}
			else
				pp.set(ptoff, p_next);

		}
	}
}
void gridUpdate::step_solution_imp(GU_Detail *gdp, Data_struct &ds, const Wilson &wil, const Material &material)
{
	GA_RWHandleV3D accp(gdp->findAttribute(GA_ATTRIB_POINT, "acc"));
	GA_RWHandleV3D vp(gdp->findAttribute(GA_ATTRIB_POINT, "v"));
	GA_RWHandleV3D pp(gdp->findAttribute(GA_ATTRIB_POINT, "P1"));
	GA_RWHandleV3D pp0(gdp->findAttribute(GA_ATTRIB_POINT, "P0"));
	GA_RWHandleV3D disp(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));
	GA_RWHandleV3 winkler_test(gdp->findAttribute(GA_ATTRIB_POINT, "wt"));
	GA_RWHandleV3D sft(gdp->findAttribute(GA_ATTRIB_POINT, "surface_force_test"));

	unsigned int point_count = gdp->getNumPoints();

	if (pp0.isValid() && pp.isValid() && vp.isValid() && accp.isValid())
	{
		unsigned int i = 0;
		
		for (GA_Iterator pointit(gdp->getPointRange()); !pointit.atEnd(); ++pointit)
		{
			const GA_Offset ptoff = *pointit;
			const GA_Index ptindex = pointit.getIndex();

			rotate_degree_of_freedom(gdp, ds, ptoff, ptindex, ds.X, false);

			pp0.set(ptoff, pp(ptoff));

			if (sft.isValid())
			{
				UT_Vector3 fds(ds.global_force(ptindex), ds.global_force(ptindex + point_count) , ds.global_force(ptindex + 2*point_count));
				sft.set(ptoff, fds);
			}
			fpreal u0 = ds.X(ptindex);
			fpreal u1 = ds.X(ptindex + point_count);
			fpreal u2 = ds.X(ptindex + 2 * point_count);
			UT_Vector3D u(u0, u1, u2);

			UT_Vector3D old_v = vp.get(ptoff);
			UT_Vector3D old_acc = accp.get(ptoff);
			UT_Vector3D old_u = pp.get(ptoff);

			UT_Vector3D new_acc = (wil.a3*u + wil.a4*old_v + wil.a5*old_acc);
			UT_Vector3D new_v = old_v + wil.a6*(new_acc + old_acc);
			
			UT_Vector3D new_u = old_u + material.dt*old_v + wil.a7*(new_acc + 2 * old_acc);
			
			if (disp.isValid())
				disp.add(ptoff, new_u - old_u );
			
			pp.set(ptoff, new_u);
			if (!wil.use_static)
			{
				vp.set(ptoff, new_v);
				accp.set(ptoff, new_acc);
			}
			
			if (winkler_test.isValid())
			{
				fpreal ps0 = ds.fixed_global_force(ptindex);
				fpreal ps1 = ds.fixed_global_force(ptindex + point_count);
				fpreal ps2 = ds.fixed_global_force(ptindex + 2 * point_count);
				UT_Vector3 pre(ps0,ps1,ps2);
				winkler_test.add(ptoff, -pre);
			}
		
			i += 1;
		}
	}
}

void gridUpdate::assemble_global_laplacian(GU_Detail *gdp, const GA_Offset &primoff, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &wil)
{
	GA_ROHandleV3D lp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P0"));
	GA_ROHandleV3D pp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P1"));
	GA_ROHandleV3D p_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
	GA_ROHandleD smd(gdp->findAttribute(GA_ATTRIB_POINT, "surface_mass_density"));
	GA_ROHandleV3D displacement_p(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));
	GA_ROHandleD disp(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
	GA_RWHandleD gvf(gdp->findAttribute(GA_ATTRIB_POINT, "gravity_flux"));
	GA_ROHandleD fd(gdp->findAttribute(GA_ATTRIB_POINT, "field_divergence"));
	GA_ROHandleD element_density(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "density"));
	const double FOURPIG = 4.0*3.14159265359*6.67408*pow(10, -11);

	if (pp_handle.isValid() && lp_handle.isValid())
	{
		unsigned int point_range = gdp->getPointRange().getEntries();

		double prim_surf_area = 0;
		double a = ds.determinant / 120.0;
		double b = ds.determinant * 1 / 6.0;
		double c = ds.surf_area / 12.0;

		GA_Primitive *prim = gdp->getPrimitive(primoff);
		UT_Vector3D force_vector(0.0);

		ds.BT.postMult(ds.N, ds.pre_stress);
		UT_MatrixD pre_stress_t(0, 11, 0, 11);

		ds.pre_stress.transpose(pre_stress_t);
		ds.pre_stress.addScaledMatrix(pre_stress_t, 1);
		
		UT_Vector vector_div_result(0, 11);
		vector_div_result.zero();
	
		ds.pre_stress.postMult(ds.divergence_vector, vector_div_result);
		
		if (!ds.stashed)
			ds.non_linear_element_3D = ds.Bs;

		double grv(0.0);
		unsigned int i_1 = 0;
		unsigned int q = 0;

		for (GA_Iterator pointit(prim->getPointRange()); !pointit.atEnd(); ++pointit)
		{
			GA_Offset ptoff = *pointit;
			GA_Index index = pointit.getIndex();
			bool is_surface = false;
			if (ds.bound(0) + ds.bound(1) + ds.bound(2) + ds.bound(3) == 3)
			{
				is_surface = true;
			}
			fpreal A = 1;
			if (ds.bound(i_1))
				A = 1;

			double surface_mass_density = 0;
			if (smd.isValid() && is_surface)
				surface_mass_density = ds.bound(i_1)*(smd(ds.tetra_points_offset(i_1)) + ds.bound(0)*smd(ds.tetra_points_offset(0)) + ds.bound(1)*smd(ds.tetra_points_offset(1)) + ds.bound(2)*smd(ds.tetra_points_offset(2)) + ds.bound(3)*smd(ds.tetra_points_offset(3)));
			
			if (!ds.stashed)
				surface_mass_density *= 0.5;
			const double boundary_conditions = disp(ptoff);

			if (disp.isValid())
			{
				if (boundary_conditions == 0 )
				{		
					ds.global_force(ds.tetra_points_offset(i_1)) -= surface_mass_density*FOURPIG*c;
					ds.global_force(ds.tetra_points_offset(i_1)) -= A*FOURPIG*element_density(primoff)*abs(b)*(vector_div_result(3*i_1) + vector_div_result(3 * i_1+1) + vector_div_result(3 * i_1+2));
					grv =  FOURPIG*element_density(primoff)*abs(b)*(vector_div_result(3 * i_1) + vector_div_result(3 * i_1 + 1) + vector_div_result(3 * i_1 + 2));
				}
			}
			if (!ds.stashed)
			{
				for (unsigned int w = 0; w < 4; w++)
				{
					inject_element_matrix_laplacian(gdp, w, q, ds, sparse, wil);
				}
			}

			if (gvf.isValid())
				gvf.add(ds.tetra_points_offset(i_1), grv);

			i_1 += 1;
			q += 1;
		}
	}
}
void gridUpdate::assemble_B_laplacian(Data_struct &ds)
{
	ds.Bs.zero();
	ds.BTs.zero();
	UT_MatrixD dx(0,0,0,3);
	UT_MatrixD dx_t(0,3,0,0);

	UT_MatrixD dy(0, 0, 0, 3);
	UT_MatrixD dy_t(0, 3, 0, 0);

	UT_MatrixD dz(0, 0, 0, 3);
	UT_MatrixD dz_t(0, 3, 0, 0);

	UT_MatrixD element_x(0, 3, 0, 3);
	UT_MatrixD element_y(0, 3, 0, 3);
	UT_MatrixD element_z(0, 3, 0, 3);

	dx.zero();
	dx_t.zero();
	dy.zero();
	dy_t.zero();
	dz.zero();
	dz_t.zero();

	element_x.zero();
	element_y.zero();
	element_z.zero();

	dx(0, 0) = ds.inverse_jacobian(0, 1);
	dx(0, 1) = ds.inverse_jacobian(1, 1);
	dx(0, 2) = ds.inverse_jacobian(2, 1);
	dx(0, 3) = ds.inverse_jacobian(3, 1);

	dy(0, 0) = ds.inverse_jacobian(0, 2);
	dy(0, 1) = ds.inverse_jacobian(1, 2);
	dy(0, 2) = ds.inverse_jacobian(2, 2);
	dy(0, 3) = ds.inverse_jacobian(3, 2);

	dz(0, 0) = ds.inverse_jacobian(0, 3);
	dz(0, 1) = ds.inverse_jacobian(1, 3);
	dz(0, 2) = ds.inverse_jacobian(2, 3);
	dz(0, 3) = ds.inverse_jacobian(3, 3);

	dx.transpose(dx_t);
	dy.transpose(dy_t);
	dz.transpose(dz_t);
	dx_t.postMult(dx, element_x);
	dy_t.postMult(dy, element_y);
	dz_t.postMult(dz, element_z);
	
	ds.Bs.addScaledMatrix(element_x,1);
	ds.Bs.addScaledMatrix(element_y,1);
	ds.Bs.addScaledMatrix(element_z,1);
	
};

void gridUpdate::inject_element_matrix_laplacian(const GU_Detail *gdp, unsigned int i, unsigned int j, Data_struct &ds, SIM_Data_Sparse *sparse, const Wilson &wil)
{
	GA_ROHandleV3 pin_p(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
	unsigned int point_count = gdp->getPointRange().getEntries();
	fpreal value = 0;
	value = ds.determinant * 1 / 6.0*(ds.non_linear_element_3D(i, j));

	unsigned int x = ds.tetra_points(i);
	unsigned int y = ds.tetra_points(j);

	if (pin_p.isValid())
	{ 
		UT_Vector3 boundsx = pin_p.get(x);
		UT_Vector3 boundsy = pin_p.get(y);

		if (boundsx(0))
		{
			if (x != y)
			{
				ds.fixed_global_force(x) -= value;
			}
			value = 0;
		}

		if (boundsy(0) == 1)
		{
			value = 0;
		}
	}
	sparse->stiffness_matrix->addToElement(x, y, value);
}
void gridUpdate::set_diagonal(Data_struct &ds, SIM_Data_Sparse *sparse, const GU_Detail *gdp)
{
	GA_ROHandleV3 disp(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
	unsigned int pointCount = gdp->getNumPoints();

	for (GA_Iterator point_iter(gdp->getPointRange()); !point_iter.atEnd(); ++point_iter)
	{
		const GA_Offset offset = *point_iter;
		const GA_Index index = point_iter.getIndex();

		if (disp.isValid())
		{
			UT_Vector3 di = disp(offset);
			if (di(0) == 1)
			{
				sparse->stiffness_matrix->addToElement(index, index, 1);
				ds.fixed_global_force(index) = 0;
			}
		}
	}
}

void gridUpdate::laplacian_solution(GU_Detail *gdp, Data_struct &ds)
{
	GA_RWHandleD pp0(gdp->findAttribute(GA_ATTRIB_POINT, "potential"));
	unsigned int point_count = gdp->getNumPoints();

	if (pp0.isValid())
	{
		unsigned int i = 0;
		for (GA_Iterator pointit(gdp->getPointRange()); !pointit.atEnd(); ++pointit)
		{
			GA_Offset ptoff = *pointit;
			GA_Index ptindex = pointit.getIndex();
			fpreal phi = ds.X(ptindex);
			pp0.set(ptoff, phi);	
		}
	}
}

void gridUpdate::calculate_and_apply_winkler(const GA_Detail *gdp, double &surface_force, const Data_struct &ds, const unsigned int &i)
{
	GA_ROHandleF density_contrast_p(gdp->findAttribute(GA_ATTRIB_POINT, "density_contrast"));
	fpreal point_range = gdp->getPointRange().getEntries();
	GA_ROHandleV3D disp(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));
	GA_ROHandleI winklerp(gdp->findAttribute(GA_ATTRIB_POINT, "winkler"));
	GA_ROHandleD potential_p(gdp->findAttribute(GA_ATTRIB_POINT, "potential"));
	GA_ROHandleV3D pp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));

	int plane_sum = ds.winkler(0) + ds.winkler(1) + ds.winkler(2) + ds.winkler(3);
	bool is_surface = false;
	if (ds.bound(0) + ds.bound(1) + ds.bound(2) + ds.bound(3) == 3)
		is_surface = true;
	bool is_core = false;
	if (ds.core(0) + ds.core(1) + ds.core(2) + ds.core(3) == 3)
		is_core = true;
	else
		is_core = false;

	if ( disp.isValid() && winklerp.isValid() && ds.core(i) && ds.core(0) + ds.core(1) + ds.core(2) + ds.core(3) == 3)
	{
		GA_Offset idx = ds.tetra_points_offset(i);
		GA_Offset id0 = ds.tetra_points_offset(0);
		GA_Offset id1 = ds.tetra_points_offset(1);
		GA_Offset id2 = ds.tetra_points_offset(2);
		GA_Offset id3 = ds.tetra_points_offset(3);

		double density_contrast = density_contrast_p(idx);
	
		UT_Vector3D p = pp_handle(idx);
		UT_Vector3D p_normal = p / p.length();
		UT_Vector3D displacement = disp(idx);
		double vertical_displacement = dot(p_normal, displacement);
		double b_force_i = ds.core(i)*vertical_displacement* 10750;
	

		density_contrast = density_contrast_p(id0);
		p = pp_handle(id0);
		p_normal = p / p.length();
		displacement = disp(id0);
		vertical_displacement = dot(p_normal, displacement);
		double b_force_0 = ds.core(0)*vertical_displacement* 10750;


		density_contrast = density_contrast_p(id1);
		p = pp_handle(id1);
		p_normal = p / p.length();
		displacement = disp(id1);
		vertical_displacement = dot(p_normal, displacement);
		double b_force_1 = ds.core(1)*vertical_displacement* 10750;


		density_contrast = density_contrast_p(id2);
		p = pp_handle(id2);
		p_normal = p / p.length();
		displacement = disp(id2);
		vertical_displacement = dot(p_normal, displacement);
		double b_force_2 = ds.core(2)*vertical_displacement* 10750;
		

		density_contrast = density_contrast_p(id3);
		p = pp_handle(id3);
		p_normal = p / p.length();
		displacement = disp(id3);
		vertical_displacement = dot(p_normal, displacement);
		double b_force_3 = ds.core(3)*vertical_displacement* 10750;
	
		surface_force += 0*ds.core(i)*(b_force_i + b_force_0 + b_force_1 + b_force_2 + b_force_3);
	}
}

void gridUpdate::calculate_and_apply_winkler(const GA_Detail *gdp, UT_Vector3D &surface_force, const Data_struct &ds, const unsigned int &i)
{
	fpreal point_range = gdp->getPointRange().getEntries();
	GA_ROHandleV3D disp(gdp->findAttribute(GA_ATTRIB_POINT, "displacement"));
	GA_ROHandleI winklerp(gdp->findAttribute(GA_ATTRIB_POINT, "winkler"));
	GA_ROHandleV3D pp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
	GA_ROHandleD gvp_w(gdp->findAttribute(GA_ATTRIB_POINT, "potential"));
	GA_ROHandleF density_contrast_p(gdp->findAttribute(GA_ATTRIB_POINT, "density_contrast"));

	if (density_contrast_p.isValid() && gvp_w.isValid() &&  ds.core(i) && ds.core(0) + ds.core(1) + ds.core(2) + ds.core(3) == 3)
	{
		GA_Offset idx = ds.tetra_points_offset(i);
		GA_Offset id0 = ds.tetra_points_offset(0);
		GA_Offset id1 = ds.tetra_points_offset(1);
		GA_Offset id2 = ds.tetra_points_offset(2);
		GA_Offset id3 = ds.tetra_points_offset(3);

		UT_Vector3D p = pp_handle(idx);
		UT_Vector3D p_normal = p / p.length();
		UT_Vector3D displacement = disp(idx);
		UT_Vector3D vertical_displacement = dot(p_normal, displacement)*p_normal;
		double g = -10.457;

		UT_Vector3D b_force_i = 10750*ds.core(i)*(vertical_displacement*g + 0*p_normal*gvp_w(idx));

		p = pp_handle(id0);
		p_normal = p / p.length();
		displacement = disp(id0);
		vertical_displacement = dot(p_normal, displacement)*p_normal;
		UT_Vector3D b_force_0 = 10750 * ds.core(0)*(vertical_displacement*g + 0*p_normal*gvp_w(id0));

		p = pp_handle(id1);
		p_normal = p / p.length();
		displacement = disp(id1);
		vertical_displacement = dot(p_normal, displacement)*p_normal;
		UT_Vector3D b_force_1 = 10750 * ds.core(1)*(vertical_displacement*g + 0*p_normal*gvp_w(id1));

		p = pp_handle(id2);
		p_normal = p / p.length();
		displacement = disp(id2);
		vertical_displacement = dot(p_normal, displacement)*p_normal;
		UT_Vector3D b_force_2 = 10750 * ds.core(2)*(vertical_displacement*g + 0*p_normal*gvp_w(id2));

		p = pp_handle(id3);
		p_normal = p / p.length();
		displacement = disp(id3);
		vertical_displacement = dot(p_normal, displacement)*p_normal;
		UT_Vector3D b_force_3 = 10750 * ds.core(3)*(vertical_displacement *g + 0*p_normal*gvp_w(id3));
			
		surface_force += ds.core(i)*(b_force_i + b_force_0 + b_force_1 + b_force_2 + b_force_3);
	}
}

void gridUpdate::rotate_degree_of_freedom(GU_Detail *gdp, GA_Offset primoff, Data_struct &ds, const unsigned int &i, const bool stashed)
{
	GA_ROHandleV3D pin(gdp->findAttribute(GA_ATTRIB_POINT, "pintoanimation"));
	GA_ROHandleD element_density(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "density"));
	GA_ROHandleD element_g(gdp->findAttribute(GA_ATTRIB_PRIMITIVE, "g0"));

	UT_Matrix rotation_matrix;
	rotation_matrix.init(0, 11, 0, 11);
	rotation_matrix.makeIdentity();

	UT_Matrix diagonal_filter(0,11,0,11);
	diagonal_filter.zero();

	UT_Matrix constraint_matrix;
	constraint_matrix.init(0, 11, 0, 11);
	constraint_matrix.zero();

	UT_Vector3D e0(1, 0, 0);
	UT_Vector3D e1(0, 1, 0);
	UT_Vector3D e2(0, 0, 1);

	UT_Vector3D r0 = gdp->getPos3(ds.tetra_points_offset(i));
	r0.normalize();
	UT_Vector3D ra(0, 0, 0);

	for (unsigned int k = 0; k < 4; k++)
	{
		GA_Offset k_off = gdp->pointOffset(k);
		UT_Vector3D candidate_pos = gdp->getPos3(5 + k_off);
		candidate_pos.normalize();
		if (dot(candidate_pos, r0) < 0.9999999)
		{
			ra = candidate_pos;
			break;
		}
	}

	UT_Vector3D r1 = cross(r0, ra);
	r1.normalize();
	UT_Vector3D r2 = cross(r0, r1);
	r2.normalize();

	if (dot(cross(r1, r2), r0) < 0)
		r2 *= -1;
	
	for (unsigned int q = 0; q < 3; q++)
	{
		for (unsigned int t = 0; t < 3; t++)
		{
			UT_Vector3D target_vector_e(0.0);
			UT_Vector3D target_vector_r(0.0);
			if (q == 0)
				target_vector_e = e0;
			if (q == 1)
				target_vector_e = e1;
			if (q == 2)
				target_vector_e = e2;
			if (t == 0)
				target_vector_r = r0;
			if (t == 1)
				target_vector_r = r1;
			if (t == 2)
				target_vector_r = r2;

			double value = dot(target_vector_e, target_vector_r);
			rotation_matrix(3 * i + q, 3 * i + t) = value;

		}
	}

	if (pin.isValid())
	{
		if (pin(ds.tetra_points_offset(i)).length() != 0)
		{
			constraint_matrix(3 * i + 1, 3 * i + 1) = 100000000.0;
			constraint_matrix(3 * i + 2, 3 * i + 2) = 100000000.0;
		}
	}
	UT_MatrixD K(0, 11, 0, 11);
	K.zero();
	K = ds.non_linear_element_3D;
	UT_MatrixD KR(0, 11, 0, 11);
	KR.zero();
	UT_MatrixD RKR(0, 11, 0, 11);
	RKR.zero();
	UT_MatrixD rotation_matrix_t(0, 11, 0, 11);

	rotation_matrix.transpose(rotation_matrix_t);
	K.postMult(rotation_matrix, KR);
	rotation_matrix_t.postMult(KR, RKR);
	ds.non_linear_element_3D = RKR;
	ds.non_linear_element_3D.addScaledMatrix(constraint_matrix, 100);
}

void gridUpdate::rotate_degree_of_freedom(const GA_Detail *gdp, Data_struct &ds, const GA_Offset &ptoff, const GA_Index &idx, UT_VectorD &target_vector, bool transpose)
{
	const unsigned int npt = gdp->getNumPoints();
	UT_Matrix3D rotation_matrix(0.0);

	UT_Vector3D e0(1, 0, 0);
	UT_Vector3D e1(0, 1, 0);
	UT_Vector3D e2(0, 0, 1);

	UT_Vector3D r0 = gdp->getPos3(ptoff);
	r0.normalize();
	UT_Vector3D ra(0, 0, 0);

	for (unsigned int k = 0; k < 4; k++)
	{
		GA_Offset k_off = gdp->pointOffset(k);
		UT_Vector3D candidate_pos = gdp->getPos3(5 + k_off);
		candidate_pos.normalize();
		if (dot(candidate_pos, r0) < 0.9999999)
		{
			ra = candidate_pos;
			break;
		}
	}

	UT_Vector3D r1 = cross(r0, ra);
	r1.normalize();
	UT_Vector3D r2 = cross(r0, r1);
	r2.normalize();

	if (dot(cross(r1, r2), r0) < 0)
		r2 *= -1;

	for (unsigned int q = 0; q < 3; q++)
	{
		for (unsigned int t = 0; t < 3; t++)
		{
			UT_Vector3D target_vector_e;
			UT_Vector3D target_vector_r;
			if (q == 0)
				target_vector_e = e0;
			if (q == 1)
				target_vector_e = e1;
			if (q == 2)
				target_vector_e = e2;

			if (t == 0)
				target_vector_r = r0;
			if (t == 1)
				target_vector_r = r1;
			if (t == 2)
				target_vector_r = r2;

			double value = dot(target_vector_e, target_vector_r);
			if(!transpose)
				rotation_matrix(t, q) = value;
			else
				rotation_matrix(q, t) = value;
		}
	}

	UT_Vector3D force_values(0,0,0);
	force_values(0) = target_vector(idx);
	force_values(1) = target_vector(npt + idx);
	force_values(2) = target_vector(2 * npt + idx);
	
	UT_Vector3D rotate_force_values = force_values*rotation_matrix;
			
	target_vector(idx) = rotate_force_values(0);
	target_vector(npt + idx) = rotate_force_values(1);
	target_vector(2 * npt + idx) = rotate_force_values(2);
	
}


void gridUpdate::add_winkler(const GA_Detail *gdp, double c, const UT_Vector3D &normal, Data_struct &ds, const unsigned int &i)
{
	GA_ROHandleF density_contrast_p(gdp->findAttribute(GA_ATTRIB_POINT, "density_contrast"));
	GA_ROHandleF gt_p(gdp->findAttribute(GA_ATTRIB_POINT, "g0"));
	fpreal point_range = gdp->getPointRange().getEntries();
	GA_ROHandleI winklerp(gdp->findAttribute(GA_ATTRIB_POINT, "winkler"));
	GA_ROHandleV3D pp_handle(gdp->findAttribute(GA_ATTRIB_POINT, "P"));
	GA_ROHandleD gvp_w(gdp->findAttribute(GA_ATTRIB_POINT, "potential"));

	if (gvp_w.isValid() && winklerp.isValid() && ds.winkler(0) + ds.winkler(1) + ds.winkler(2) + ds.winkler(3) == 3 && density_contrast_p.isValid() && gt_p.isValid())
	{
		GA_Offset idx = ds.tetra_points_offset(i);
		GA_Offset id0 = ds.tetra_points_offset(0);
		GA_Offset id1 = ds.tetra_points_offset(1);
		GA_Offset id2 = ds.tetra_points_offset(2);
		GA_Offset id3 = ds.tetra_points_offset(3);

		UT_Vector3D pos_i = gdp->getPos3(idx);
		UT_Vector3D n_p = pos_i / pos_i.length();
		double scaled_c = c*abs(dot(n_p, normal));

		double density_contrast = density_contrast_p(idx);
		double g = -(gt_p(idx));
		double b_force_i = density_contrast_p(idx)*(g*ds.winkler(i) );


		density_contrast = density_contrast_p(id0);
		g = -(gt_p(id0));
		double b_force_0 = density_contrast_p(id0) *(g*ds.winkler(0));


		density_contrast = density_contrast_p(id1);
		g = -(gt_p(id1));
		double b_force_1 = density_contrast_p(id1)*(g*ds.winkler(1));


		density_contrast = density_contrast_p(id2);
		g = -(gt_p(id2));
		double b_force_2 = density_contrast_p(id2)*(g*ds.winkler(2));


		density_contrast = density_contrast_p(id3);
		g = -(gt_p(id3));
		double b_force_3 = density_contrast_p(id3) *(g*ds.winkler(3));

		double winkler_addon = (6.0*scaled_c / ds.determinant)*ds.winkler(i)*(b_force_i + b_force_0 + b_force_1 + b_force_2 + b_force_3);
		ds.non_linear_element_3D(3 * i, 3 * i) -= winkler_addon;
	}
}
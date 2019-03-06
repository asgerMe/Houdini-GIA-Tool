#include "SIM_FEMSolver.h"
#include "Solve_System.h"
#include <UT/UT_DSOVersion.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_SingleSolver.h>
#include <SIM/SIM_GeometryCopy.h>
#include <UT/UT_SparseMatrix.h>
#include "Preconditioner.h"
#include "IterationTester.h"
#include <GU/GU_DetailHandle.h>
#include <UT/UT_Matrix.h>
#include <UT/UT_MatrixSolver.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_ObjectArray.h>
#include <SIM/SIM_Engine.h>
#include "GETFEMMatrices.h"
#include "Placeholder.h"
#include <GU/GU_SDF.h>
#include <PRM/PRM_ChoiceList.h>
#include <GAS/GAS_SubSolver.h>


using namespace HDK_Sample;
void
initializeSIM(void *)
{
	IMPLEMENT_DATAFACTORY(SIM_FEMSolver);
	IMPLEMENT_DATAFACTORY(SIM_Data_Sparse);
}
SIM_FEMSolver::SIM_FEMSolver(const SIM_DataFactory *factory) : BaseClass(factory)
{
}
SIM_FEMSolver::~SIM_FEMSolver()
{
}

const SIM_DopDescription *
SIM_FEMSolver::getFEMSolverDopDescription()
{
	static PRM_Name theGeometryName(GAS_NAME_GEOMETRY, "Geometry");
	static PRM_Name      theConName(SIM_NAME_CONSTRAINTPASSES, "Sub Steps");
	static PRM_Name         sopOrdinalName(SIM_NAME_INTEGRATION_MODE, "Integration Mode");
	static PRM_Name			useSavedData(SIM_NAME_USESAVED, "Use Saved Data structs");
	static PRM_Name         poissonName(SIM_NAME_POISSON, "Poisson Ratio");
	static PRM_Name         densityName(SIM_NAME_DENSITY, "Density");
	static PRM_Name         viscosityName(SIM_NAME_VISCOSITY, "Viscosity");
	static PRM_Name         youngName(SIM_NAME_YOUNG, "Young's Modulus");
	static PRM_Name			scaleTimeName(SIM_NAME_SCALETIME, "Scale Time");
	
	static PRM_Default youngDefault(100);
	static PRM_Default poissonDefault(0.25);
	static PRM_Default densityDefault(1);
	static PRM_Default viscosityDefault(0);
	static PRM_Default substepDefault(1);

	static PRM_Default winklerDefault(0);

	static PRM_Range viscosityRange(PRM_RANGE_RESTRICTED, 0);
	static PRM_Range densityRange(PRM_RANGE_RESTRICTED, 0.01);
	static PRM_Range poissonRange(PRM_RANGE_RESTRICTED, 0.001, PRM_RANGE_RESTRICTED, 0.499);
	static PRM_Range youngsRange(PRM_RANGE_RESTRICTED, 0.01);


	static PRM_Name         sopOrdChoices[] =
	{
		PRM_Name("choice2", "Implicit"),
		PRM_Name("choice3", "Explicit"),
		PRM_Name("choice4", "Static"),
		PRM_Name("choice5", "Solve Laplacian"),
		PRM_Name(0)
	};
	static PRM_ChoiceList   sopOrdinalMenu(PRM_CHOICELIST_SINGLE, sopOrdChoices);

	static PRM_Template          theTemplates[] = {
		PRM_Template(PRM_STRING, 1, &theGeometryName),
		PRM_Template(PRM_TOGGLE, 1, &useSavedData, PRMzeroDefaults),
		PRM_Template(PRM_FLT, 1, &scaleTimeName,  PRMoneDefaults, 0, &viscosityRange),
		PRM_Template(PRM_INT,1,&theConName, &substepDefault),
		PRM_Template(PRM_ORD, 1, &sopOrdinalName, 0, &sopOrdinalMenu),
		PRM_Template(PRM_FLT, 1, &youngName,  &youngDefault, 0, &youngsRange),
		PRM_Template(PRM_FLT, 1, &poissonName, &poissonDefault, 0, &poissonRange),
		PRM_Template(PRM_FLT, 1, &densityName, &densityDefault, 0, &densityRange),
		PRM_Template(PRM_FLT, 1, &viscosityName, &viscosityDefault, 0, &viscosityRange),
		PRM_Template()
	};

	static SIM_DopDescription    theDopDescription(true,
		"gas_visco_elastic_fem",
		"GAS VEFEM_Solver",
		"Solver",
		classname(),
		theTemplates);
	setGasDescription(theDopDescription);
	return &theDopDescription;
}

bool SIM_FEMSolver::solveGasSubclass(
	SIM_Engine& engine, SIM_Object* object,
	SIM_Time time,
	SIM_Time time_step)
{
	SIM_GeometryCopy *state = getGeometryCopy(object, GAS_NAME_GEOMETRY);
	
	if (!state)
	{
		addError(object, SIM_MISSINGDATA,
			"Geometry", UT_ERROR_MESSAGE);
		
		return false;
	}

	if (state)
	{	
		unsigned int sub_steps = getConstraintPasses();
		GU_DetailHandleAutoWriteLock lock(state->lockGeometry());
		GU_Detail *gdp = lock.getGdp();
		SIM_DataArray data_array;
		if (gdp->getNumPoints() == 0)
			return false;
		Material parms;
		parms.density = getdensity();
		parms.viscosity = getviscosity();
		parms.youngs_modulus = getyoung();
		parms.poisson_ratio = getpoisson();
		
		parms.dt = time_step;
		parms.dt *= getscaletime();
		UT_String data_name;
		float  theta = 1.7;
	
		int __explicit__ = 0;
		if (getintMode() == 0)
		{
			__explicit__ = 0;
			data_name = "SparseData";
		}
		if (getintMode() == 1)
			__explicit__ = 1;
		if (getintMode() == 2)
		{
			__explicit__ = 0;
			theta = 1;
			data_name = "SparseData";
		}
		if (getintMode() == 3)
		{
			__explicit__ = -2;
			data_name = "SparseData_Laplacian";
		}

		getMatchingDataByName(data_array, object, data_name);
		SIM_Data_Sparse *sparse_data;
		if (data_array.entries() != 0)
		{
			sparse_data = SIM_DATA_CAST(data_array[0], SIM_Data_Sparse);
			Wilson w(parms.dt, theta);
			SOLVE_LINEAR::SOLVE_3D(sparse_data, gdp, parms, w, sub_steps, __explicit__);
			gdp->bumpAllDataIds();
			pubHandleModification();
		}
		else
		{
			addError(object, SIM_MISSINGDATA,
				"Geometry", UT_ERROR_MESSAGE);
			return false;
		}

	}
	return state ? true : false;
}



SIM_Data_Sparse::SIM_Data_Sparse(const SIM_DataFactory *factory) : BaseClass(factory), stiffness_matrix(0), fixed_global(0)
{

}
SIM_Data_Sparse::~SIM_Data_Sparse()
{
	freeArray();
}

const SIM_DopDescription *
SIM_Data_Sparse::getSparseDopDescription()
{
	static PRM_Template          theTemplates[] = { PRM_Template() };

	static SIM_DopDescription    theDopDescription(true,
		"sparse_dataholder",
		"Sparse Data Container",
		"SparseData",
		classname(),
		theTemplates);
	return &theDopDescription;
}
void
SIM_Data_Sparse::initializeSubclass()
{
	BaseClass::initializeSubclass();
	freeArray();
	stiffness_matrix = 0;
	fixed_global = 0;
}

void
SIM_Data_Sparse::makeEqualSubclass(const SIM_Data *source)
{
	BaseClass::makeEqualSubclass(source);
	
	const SIM_Data_Sparse     *srcvox;
	srcvox = SIM_DATA_CASTCONST(source, SIM_Data_Sparse);
	if (srcvox)
	{
		if (srcvox->stiffness_matrix)
		{
			allocateArray();
			setDivisions(srcvox->getDivisions());
			*stiffness_matrix = *srcvox->stiffness_matrix;
		}
		else
		{
			freeArray();
		}
	}
}

void SIM_Data_Sparse::freeArray() const
{
	delete stiffness_matrix;
	delete fixed_global;

	stiffness_matrix = 0;
	fixed_global = 0;
}

void SIM_Data_Sparse::allocateArray() const
{
	if (!stiffness_matrix)
		stiffness_matrix = new UT_SparseMatrix;

	if (!fixed_global)
		fixed_global = new UT_Vector;
}

void SIM_Data_Sparse::setDivisions(unsigned int size)
{
	if (!stiffness_matrix || !fixed_global)
		allocateArray();
	
	stiffness_matrix->init(size, size);

	fixed_global->init(0, size-1);
	fixed_global->zero();
}

unsigned int SIM_Data_Sparse::getDivisions() const
{
	if(!stiffness_matrix)
		allocateArray();
	return stiffness_matrix->getNumRows();
}


void
SIM_Data_Sparse::saveIOSubclass(std::ostream &os, SIM_DataThreadedIO *io) const
{

	BaseClass::saveIOSubclass(os, io);
	stiffness_matrix->save(os);
}
bool
SIM_Data_Sparse::loadIOSubclass(UT_IStream &is, SIM_DataThreadedIO *io)
{
	if (!BaseClass::loadIOSubclass(is, io))
		return false;
	stiffness_matrix->load(is);
	return true;
}

int64 SIM_Data_Sparse::getMemorySizeSubclass() const
{
	if (stiffness_matrix)
		return	stiffness_matrix->getMemoryUsage();
	else
		return 0;
}

bool SIM_Data_Sparse::is_loaded() const 
{
	if (stiffness_matrix)
	{
		if (stiffness_matrix->getMemoryUsage() > int64(1000))
			return true;
		else
			return false;
	}
	else
		return false;
}


#ifndef __SIM_FEMSOLVER_h__
#define __SIM_FEMSOLVER_h__

#include <SIM/SIM_OptionsUser.h>
#include <SIM/SIM_SingleSolver.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SIM/SIM_Utils.h>
#include "Placeholder.h"
#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

#define SIM_NAME_CONSTRAINTPASSES     "cPasses"
#define SIM_NAME_ERRORTOLERANCE "errorTol"
#define SIM_NAME_INTEGRATION_MODE "intMode"

#define SIM_NAME_YOUNG "Y"
#define SIM_NAME_POISSON "v"
#define SIM_NAME_DENSITY "density"
#define SIM_NAME_VISCOSITY "viscosity"

#define SIM_NAME_SCALETIME "scaletime"
#define GAS_NAME_PATH "path"
#define SIM_NAME_USESAVED "useSaved"



namespace HDK_Sample {
	class SIM_FEMSolver : public GAS_SubSolver
	{
	public:
		GETSET_DATA_FUNCS_I(SIM_NAME_CONSTRAINTPASSES, ConstraintPasses);
		GETSET_DATA_FUNCS_F(SIM_NAME_ERRORTOLERANCE, errorTolerance);
		GETSET_DATA_FUNCS_I(SIM_NAME_INTEGRATION_MODE, intMode);
		GETSET_DATA_FUNCS_F(SIM_NAME_SCALETIME, scaletime);
		GETSET_DATA_FUNCS_F(SIM_NAME_YOUNG, young);
		GETSET_DATA_FUNCS_F(SIM_NAME_POISSON, poisson);
		GETSET_DATA_FUNCS_F(SIM_NAME_DENSITY, density);
		GETSET_DATA_FUNCS_F(SIM_NAME_VISCOSITY, viscosity);
		GETSET_DATA_FUNCS_S(GAS_NAME_PATH, matrixpath);
		GETSET_DATA_FUNCS_I(SIM_NAME_USESAVED, usesaved);
	
		void pubHandleModification(){
			handleModification();
		}

	protected:

		explicit             SIM_FEMSolver(const SIM_DataFactory *factory);
		virtual             ~SIM_FEMSolver();

		virtual bool solveGasSubclass(
			SIM_Engine& engine, SIM_Object* object,
			SIM_Time time,
			SIM_Time time_step);

	private:
		static const SIM_DopDescription     *getFEMSolverDopDescription();
		DECLARE_STANDARD_GETCASTTOTYPE();
		DECLARE_DATAFACTORY(SIM_FEMSolver,
			GAS_SubSolver,
			"GAS VEFEM_Solver",
			getFEMSolverDopDescription());
	};
} // End HDK_Sample namespace
class SIM_Data_Sparse : public SIM_Data
{
protected:
	explicit             SIM_Data_Sparse(const SIM_DataFactory *factory);
	virtual             ~SIM_Data_Sparse();

	virtual void         initializeSubclass();
	virtual void         makeEqualSubclass(const SIM_Data *source);
	virtual bool loadIOSubclass(UT_IStream &is, SIM_DataThreadedIO *io);
	virtual void saveIOSubclass(std::ostream &os, SIM_DataThreadedIO *io) const;
	virtual int64        getMemorySizeSubclass() const;

public:
	void pubHandleModification() { handleModification(); };
	void setDivisions(unsigned int size);
	unsigned int getDivisions() const;
	bool is_loaded() const;
	mutable UT_SparseMatrix				*stiffness_matrix;
	mutable UT_SparseMatrix				*stiffness_matrix_s;
	mutable UT_SparseMatrix				*laplacian_stiffness_matrix;
	mutable UT_Vector *fixed_global;
	mutable UT_Vector *laplacian_fixed_global;

private:
	static const SIM_DopDescription     *getSparseDopDescription();
	void freeArray() const;
	void allocateArray() const;

	DECLARE_STANDARD_GETCASTTOTYPE();
	DECLARE_DATAFACTORY(SIM_Data_Sparse,
	SIM_Data,
		"Data Holder",
		getSparseDopDescription()
		);
};
#endif



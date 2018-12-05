
#ifndef SOFA_COMPONENT_FORCEFIELD_Friction_H
#define SOFA_COMPONENT_FORCEFIELD_Friction_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include "initPlugin.h"

#include <SofaMiscFem/StandardTetrahedralFEMForceField.h>
#include <SofaMiscFem/HyperelasticMaterial.h>
#include <sofa/core/behavior/ForceField.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/MatSym.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaBaseTopology/TopologyData.h>
#include <string>
#include <map>

#include <sofa/component/component.h>
#include <string>
#include <iostream>
#include <math.h>

#include <Eigen/src/Core/Matrix.h>

namespace sofa
{

namespace component
{

namespace forcefield
{
using namespace std;
using namespace sofa::defaulttype;
using namespace sofa::component::topology;


//***************** Tetrahedron FEM code for several elastic models: StandardTetrahedralFEMForceField*******************************************************************
//********************************** Based on classical discretization : Fi=-Bi^T S V and Kij=Bi^T N Bj +Di^T S Dj **********************************************
//***************************************** where Bi is the strain displacement (6*3 matrix), S SPK tensor N=dS/dC, Di shape vector ************************************
//**************************** Code dependant on HyperelasticMatrialFEM and inherited classes *********************************************************************

/** Compute Finite Element forces based on tetrahedral elements.
*/
template<class DataTypes>
class Friction: public StandardTetrahedralFEMForceField<DataTypes>
{
  public:
	  SOFA_CLASS(SOFA_TEMPLATE(Friction, DataTypes), SOFA_TEMPLATE(StandardTetrahedralFEMForceField, DataTypes));

    typedef StandardTetrahedralFEMForceField<DataTypes> Inherited;
    typedef typename Inherited::VecCoord VecCoord;
    typedef typename Inherited::VecDeriv VecDeriv;
    typedef typename Inherited::Coord Coord;
    typedef typename Inherited::Deriv Deriv;
    typedef typename Inherited::Real Real;
    typedef defaulttype::Mat<3,3,Real> Matrix3;
    typedef defaulttype::Mat<6,6,Real> Matrix6;
    typedef defaulttype::Mat<6,3,Real> Matrix63;
    typedef defaulttype::MatSym<3,Real> MatrixSym;

    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv; 
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord; 

    typedef helper::vector<Real> SetParameterArray;
    typedef helper::vector<Coord> SetAnisotropyDirectionArray;

    typedef core::topology::BaseMeshTopology::index_type Index;
    typedef core::topology::BaseMeshTopology::Tetra Element;
    typedef core::topology::BaseMeshTopology::SeqTetrahedra VecElement;

    typedef typename Inherited::EdgeInformation EdgeInformation;
    typedef typename Inherited::edgeInformationVector edgeInformationVector;
    typedef typename Inherited::tetrahedronRestInfoVector tetrahedronRestInfoVector;
    typedef typename Inherited::TetrahedronRestInformation TetrahedronRestInformation;

public :


protected:



public:

   Friction();	
	virtual ~Friction();
	
	//enum pointCollisionState {  COLLIDE_COLLIDE_ANCHORCHANGE_FWD, COLLIDE_COLLIDE_ANCHORCHANGE_BWD, COLLIDE_COLLIDE_NOTANCHORCHANGE_FWD, COLLIDE_COLLIDE_NOTANCHORCHANGE_BWD, NOTCOLLIDE_COLLIDE, COLLIDE_NOTCOLLIDE, NOTCOLLIDE_NOTCOLLIDE};
	enum pointCollisionState {  COLLIDE_COLLIDE_ANCHORCHANGE_1, COLLIDE_COLLIDE_ANCHORCHANGE_2, COLLIDE_COLLIDE_ANCHORCHANGE_12, COLLIDE_COLLIDE_NOTANCHORCHANGE, NOTCOLLIDE_COLLIDE, COLLIDE_NOTCOLLIDE, NOTCOLLIDE_NOTCOLLIDE};
	enum frictionElection {FORWARD, BACKWARD};
    struct nodeProperties{
		std::vector<double> m_anchor;
		double m_penalty;
		pointCollisionState m_collisionState;
		frictionElection m_frictionElection;
		frictionElection m_frictionElection2;
		std::vector<double> m_anchorDifference;
		std::vector<double> m_anchorDifferenceUnit;
		double m_differenceNorm;
		double m_differenceNorm2;

		double m_currentMu;
		double m_currentScalarSum;
		double m_currentScalarMuSum;

		std::vector<double> m_directionVector;
		double m_directionNorm;
	};

	void init();
	void reset();
	void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx);
	void addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv&  d_v );
	void addKToMatrix(const sofa::core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix);
	void draw(const core::visual::VisualParams* vparams);

	std::vector<std::string> &split(std::string &s, char delim, std::vector<std::string> &elems);

	std::vector<std::string> split(std::string &s, char delim);

	virtual bool load(const char* filename);

	
	
	bool m_firstStep;
	Data< VecCoord > m_mu;
	Data<Real> m_floorK;
	Data< VecCoord > m_springK;
	Data<Real> m_floorHeight;
	std::vector<nodeProperties> m_nodes;
	Data<bool> f_writeTraces;
	Data<bool> f_writeTracesD;
	sofa::core::objectmodel::DataFileName f_fileFriction;
	std::map<int,int> m_frictionRelations;

	std::vector<unsigned int> myList;
	std::vector<double> frogSum;
};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
#endif
#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_SIELEGANSPLUGIN)

#ifndef SOFA_FLOAT
extern template class Friction<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class Friction<Vec3fTypes>;
#endif

#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_SIELEGANSPLUGIN)


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_MySTANDARDTETRAHEDRALFEMFORCEFIELD_H

#ifndef SOFA_COMPONENT_FORCEFIELD_StandardTetrahedralFEMForceFieldTwoMaterials_H
#define SOFA_COMPONENT_FORCEFIELD_StandardTetrahedralFEMForceFieldTwoMaterials_H

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

#include <LastFriction.h>
#include <sofa/simulation/common/Node.h>

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
class StandardTetrahedralFEMForceFieldTwoMaterials: public StandardTetrahedralFEMForceField<DataTypes>
{
  public:
	  SOFA_CLASS(SOFA_TEMPLATE(StandardTetrahedralFEMForceFieldTwoMaterials, DataTypes), SOFA_TEMPLATE(StandardTetrahedralFEMForceField, DataTypes));

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
	vector<vector<double> > internalForce;
	vector<vector<double> > conservForces;
	vector<vector<double> > forceCorrective;
	int m; // Num internalForce rows
	int n; // Num internalForce columns
	bool firstime; // First loop
	bool conservationLinearMomentum;
	bool conservationAngularMomentum;
	double timestep;
	double counter;
	int row;

	typename sofa::component::fem::MaterialParameters<DataTypes> globalParameters2;
	Data<unsigned> f_secondMaterialStart;

protected:

	Data<SetParameterArray> f_secondParameterSet;


public:

   StandardTetrahedralFEMForceFieldTwoMaterials();	
	virtual ~StandardTetrahedralFEMForceFieldTwoMaterials();
    
	void init();
	void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx);
	void addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& /* d_v */);

	void draw(const core::visual::VisualParams* vparams);

	std::string intToString(int number);
	int loadVector(std::string filename, int i);
	vector<double> interpolate(int row);
	Eigen::Vector3d createForceCorrective(Eigen::Vector3d position, Eigen::Vector3d centerOfMass, Eigen::Vector3d lambdaForce,  Eigen::Vector3d lambdaTorque);
	Eigen::Vector3d centerMass(const VecCoord& x, int numVert);
	void multLagrange(Eigen::Vector3d &lambdaForce, Eigen::Vector3d &lambdaTorque, const VecCoord& x, Eigen::Vector3d force, Eigen::Vector3d torque, Eigen::Vector3d centerOfMass, int numVert);

	Data<bool> f_writeTraces;
	bool m_writeTraces;

	//drawing variables
	simulation::Node* m_gnode;
	std::map<int,std::vector<int>>* m_frictionRelations;
	std::vector<sofa::component::forcefield::LastFriction<Vec3dTypes>::nodeProperties>* m_nodes;
	VecDeriv auxf;

};

#ifndef SOFA_FLOAT
using sofa::defaulttype::Vec3dTypes;
#endif
#ifndef SOFA_DOUBLE
using sofa::defaulttype::Vec3fTypes;
#endif

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_SIELEGANSPLUGIN)

#ifndef SOFA_FLOAT
extern template class StandardTetrahedralFEMForceFieldTwoMaterials<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class StandardTetrahedralFEMForceFieldTwoMaterials<Vec3fTypes>;
#endif

#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_SIELEGANSPLUGIN)


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_StandardTetrahedralFEMForceFieldTwoMaterials_H
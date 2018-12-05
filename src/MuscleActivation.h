/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_MASS_MuscleActivation_H
#define SOFA_COMPONENT_MASS_MuscleActivation_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Event.h>
#include <SofaBaseTopology/TopologyData.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/RigidTypes.h>
//VERY IMPORTANT FOR GRAPHS
#include <sofa/helper/map.h>

#include <sofa/core/topology/BaseMeshTopology.h>


#include <SofaBaseTopology/EdgeSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
#include <SofaBaseTopology/QuadSetGeometryAlgorithms.h>
#include <SofaBaseTopology/BezierTetrahedronSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TetrahedronSetGeometryAlgorithms.h>
#include <SofaBaseTopology/HexahedronSetGeometryAlgorithms.h>

#ifdef SOFA_HAVE_EIGEN2
#include <SofaEigen2Solver/EigenSparseMatrix.h>
#endif

namespace sofa
{

namespace component
{
namespace topology
{
	/// forward declaration to avoid adding includes in .h
	template< class DataTypes> class EdgeSetGeometryAlgorithms;
	template< class DataTypes> class TriangleSetGeometryAlgorithms;
	template< class DataTypes> class TetrahedronSetGeometryAlgorithms;
	template< class DataTypes> class BezierTetrahedronSetGeometryAlgorithms;
	template< class DataTypes> class QuadSetGeometryAlgorithms;
	template< class DataTypes> class HexahedronSetGeometryAlgorithms;
}

namespace mass
{

using namespace sofa::component::topology;
using namespace sofa::defaulttype;

template<class DataTypes, class TMassType>
class MuscleActivationInternalData
{
};



// template<class Vec> void readVec1(Vec& vec, const char* str);
template <class DataTypes, class TMassType>
class MuscleActivation : public core::behavior::Mass<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(MuscleActivation,DataTypes,TMassType), SOFA_TEMPLATE(core::behavior::Mass,DataTypes));

    typedef core::behavior::Mass<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord                    VecCoord;
    typedef typename DataTypes::VecDeriv                    VecDeriv;
    typedef typename DataTypes::Coord                       Coord;
    typedef typename DataTypes::Deriv                       Deriv;
    typedef typename DataTypes::Real                        Real;
    typedef core::objectmodel::Data<VecCoord>               DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>               DataVecDeriv;
    typedef TMassType                                       MassType;
    typedef helper::vector<MassType>                        MassVector;
    typedef helper::vector<MassVector>                        MassVectorVector;

    // In case of non 3D template
    typedef Vec<3,Real> Vec3;
    typedef StdVectorTypes< Vec3, Vec3, Real >     GeometricalTypes ; /// assumes the geometry object type is 3D

    /// Topological enum to classify encounter meshes
    typedef enum
    {
        TOPOLOGY_UNKNOWN=0,
        TOPOLOGY_EDGESET=1,
        TOPOLOGY_TRIANGLESET=2,
        TOPOLOGY_TETRAHEDRONSET=3,
        TOPOLOGY_QUADSET=4,
        TOPOLOGY_HEXAHEDRONSET=5,
        TOPOLOGY_BEZIERTETRAHEDRONSET=6,
    } TopologyType;
	/// the way the mass should be computed on non-linear elements
	typedef enum 
	{
		EXACT_INTEGRATION=1,
		NUMERICAL_INTEGRATION=2,
		AFFINE_ELEMENT_INTEGRATION=3
	} IntegrationMethod;


    /// Mass info are stocked on vertices and edges (if lumped matrix)
    PointData<helper::vector<MassType> >  vertexMassInfo;
    EdgeData<helper::vector<MassType> >   edgeMassInfo;

    /* ---------- Specific data for Bezier Elements ------*/
    /// use this data structure to store mass for Bezier tetrahedra. 
    //// The size of the vector is nbControlPoints*(nbControlPoints+1)/2 where nbControlPoints=(degree+1)*(degree+2)*(degree+3)/2
    TetrahedronData<helper::vector<MassVector> > tetrahedronMassInfo;
    // array of Tetrahedral Bezier indices
    //sofa::helper::vector<TetrahedronBezierIndex> tbiArray;
    /* ---------- end ------*/

    /// the mass density used to compute the mass from a mesh topology and geometry
    Data< Real >         m_massDensity;

    /// to display the center of gravity of the system
    Data< bool >         showCenterOfGravity;
    Data< Real >         showAxisSize;
    /// if mass lumping should be performed (only compute mass on vertices)
    Data< bool >         lumping;
    /// if specific mass information should be outputed
    Data< bool >         printMass;
    Data<std::map < std::string, sofa::helper::vector<double> > > f_graph;
    /// the order of integration for numerical integration
    Data<size_t>	     numericalIntegrationOrder;
    /// the type of numerical integration method chosen
    Data<size_t>	     numericalIntegrationMethod;
    /// the type of integration method chosen for non linear element.
    Data<std::string>	 d_integrationMethod; 
    IntegrationMethod    integrationMethod;


	sofa::core::objectmodel::DataFileName f_forceFiles;
	sofa::core::objectmodel::DataFileName f_fileNodes;
	sofa::core::objectmodel::DataFileName f_fileActivation;
	sofa::core::objectmodel::DataFileName f_signalFile;
	sofa::core::objectmodel::DataFileName f_musclePositionFile;
	sofa::core::objectmodel::DataFileName f_fileExtremes;
	std::vector<std::vector<unsigned int>> m_activationOrder;
	std::vector<std::vector<unsigned int>> m_extremeIndices;
	std::vector<std::vector<std::vector<std::vector<unsigned int>>>> m_orderedMuscles;
	std::vector<std::vector<double>> m_activationSignals;
	std::map<std::string,std::vector<unsigned int>> m_muscleMatching;
	Real m_activationConstant1;
	Real m_activationConstant2;
	Real m_activationStop1;
	Real m_activationStop2;
	Data<Real> f_activationConstant1;
	Data<Real> f_activationConstant2;
	Data<Real> f_activationStop1;
	Data<Real> f_activationStop2;
	Real m_activationPeriod1;
	Real m_activationPeriod2;
	Data<Real> f_activationPeriod1;
	Data<Real> f_activationPeriod2;
	Data<Real> f_neuronTimestep;
	Data<int> f_minActivatedMuscles;
	int m_minActivatedMuscles;
	Data<int> f_maxActivatedMuscles;
	int m_maxActivatedMuscles;
	std::vector< std::vector< double > > m_saveSignals;
	std::vector< double > m_saveSignalsOrdered;

	Data<bool> f_correctiveForces;
	std::vector<Eigen::Vector3d> m_forceHelper;
	double m_invTotalMass;
	std::vector<Real> m_massHelper;

	Eigen::Vector3d createForceCorrective(Eigen::Vector3d position, Eigen::Vector3d centerOfMass, Eigen::Vector3d lambdaForce,  Eigen::Vector3d lambdaTorque);
	Eigen::Vector3d centerMass(const VecCoord& x, int numVert);
	void multLagrange(Eigen::Vector3d &lambdaForce, Eigen::Vector3d &lambdaTorque, const VecCoord& x, Eigen::Vector3d force, Eigen::Vector3d torque, Eigen::Vector3d centerOfMass, int numVert);



protected:

    /// The type of topology to build the mass from the topology
    TopologyType topologyType;
    Real massLumpingCoeff;
    Real savedMass;

    MuscleActivation();
    ~MuscleActivation();

    /// Internal data required for Cuda computation (copy of vertex mass for deviceRead)
    MuscleActivationInternalData<DataTypes, MassType> data;
    friend class MuscleActivationInternalData<DataTypes, MassType>;

public:

    sofa::core::topology::BaseMeshTopology* _topology;

    sofa::component::topology::EdgeSetGeometryAlgorithms<GeometricalTypes>* edgeGeo;
    sofa::component::topology::TriangleSetGeometryAlgorithms<GeometricalTypes>* triangleGeo;
    sofa::component::topology::QuadSetGeometryAlgorithms<GeometricalTypes>* quadGeo;
    sofa::component::topology::TetrahedronSetGeometryAlgorithms<GeometricalTypes>* tetraGeo;
    sofa::component::topology::HexahedronSetGeometryAlgorithms<GeometricalTypes>* hexaGeo;
    sofa::component::topology::BezierTetrahedronSetGeometryAlgorithms<GeometricalTypes>* bezierTetraGeo;

    void clear();

    virtual void reinit();
    virtual void init();

	virtual bool loadForceNodeFile(const char* filename);
	virtual bool loadActivationOrderFile(const char* filename);
	virtual bool loadSignalFile(const char* filename);
	virtual bool loadMuscleExtremes(const char* filename);
	virtual bool loadMusclePositions(const char* filename);

	void draw(const core::visual::VisualParams* vparams);

    TopologyType getMassTopologyType() const
    {
        return topologyType;
    }

    void setMassTopologyType(TopologyType t)
    {
        topologyType = t;
    }


    Real getMassDensity() const
    {
        return m_massDensity.getValue();
    }

    void setMassDensity(Real m)
    {
        m_massDensity.setValue(m);
    }

    /// Copy the vertex mass scalar (in case of CudaTypes)
    void copyVertexMass();


    // -- Mass interface
    void addMDx(const core::MechanicalParams* /* PARAMS FIRST */, DataVecDeriv& f, const DataVecDeriv& dx, double factor);

    void accFromF(const core::MechanicalParams* /* PARAMS FIRST */, DataVecDeriv& a, const DataVecDeriv& f); // This function can't be used as it use M^-1

    virtual void addForce(const core::MechanicalParams* /* PARAMS FIRST */mp, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

    double getKineticEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecDeriv& v) const;  ///< vMv/2 using dof->getV()

    double getPotentialEnergy(const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& x) const;   ///< Mgx potential in a uniform gravity field, null at origin

    virtual defaulttype::Vec6d getMomentum(const core::MechanicalParams* mparams /* PARAMS FIRST */, const DataVecCoord& x, const DataVecDeriv& v) const;  ///< (Mv,cross(x,Mv))

    void addGravityToV(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_v);

    bool isDiagonal() {return false;}

	std::vector<std::string> &split(std::string &s, char delim, std::vector<std::string> &elems);

	std::vector<std::string> split(std::string &s, char delim);

	std::vector<double> getSavedSignals(){
		return m_saveSignalsOrdered;
	}


    void setMuscleActivation(std::vector<double> values)
    {
		std::cout << "set muscles" << std::endl;
		if(values.size() != 95) {
			std::cout << "ERROR: not proper muscle values" << std::endl;
			int sz = (int) values.size();
			printf("These are the received values: (size %u)\n",sz);
			//std::cout << "These are the received values: (size = " << i <<")" << std::endl;
			//for(int i = 0; i < sz; i++){
			//	printf("%lf ", values[i]);
			//	std::cout << values.at(i) << " ";
			//}
			std::cout << std::endl;
		}
		m_activationValues = values;
     }



    /// Add Mass contribution to global Matrix assembling
    void addMToMatrix(const core::MechanicalParams *mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix);

    double getElementMass(unsigned int index) const;
    void getElementMass(unsigned int index, defaulttype::BaseMatrix *m) const;

    //void draw(const core::visual::VisualParams* vparams);

    /// Answer wether mass matrix is lumped or not
    bool isLumped() { return lumping.getValue(); }

		//amujika
	public:
		std::vector<unsigned int> m_forceNodes;
		int rowNum; // Num internalForce rows
		int colNum; // Num internalForce columns
		bool firstime; // First loop
		double timestep;
		double counter;
		int row;
		std::vector<double> m_activationValues;


protected:

    class VertexMassHandler : public topology::TopologyDataHandler<Point,MassVector>
    {
    public:
        VertexMassHandler(MuscleActivation<DataTypes,TMassType>* _m, PointData<helper::vector<TMassType> >* _data) : topology::TopologyDataHandler<Point,helper::vector<TMassType> >(_data), m(_m) {}

        /// Mass initialization Creation Functions:
        /// Vertex mass coefficient matrix creation function
        void applyCreateFunction(unsigned int pointIndex, TMassType & VertexMass,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double >&);


        ///////////////////////// Functions on Triangles //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Triangular Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new triangles
        void applyTriangleCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Triangle >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new triangles
        void applyTriangleDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesAdded* /*event*/);
        /// Callback to remove triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesRemoved* /*event*/);


        ///////////////////////// Functions on Quads //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Quad Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new quads
        void applyQuadCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Quad >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new quads
        void applyQuadDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add quads elements.
        void ApplyTopologyChange(const core::topology::QuadsAdded* /*event*/);
        /// Callback to remove quads elements.
        void ApplyTopologyChange(const core::topology::QuadsRemoved* /*event*/);


        ///////////////////////// Functions on Tetrahedron //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Tetrahedral Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new tetrahedra
        void applyTetrahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Tetrahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new tetrahedra
        void applyTetrahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraAdded* /*event*/);
        /// Callback to remove tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraRemoved* /*event*/);


        ///////////////////////// Functions on Hexahedron //////////////////////////////////////

        /// Mass coefficient Creation/Destruction functions for Hexahedral Mesh:
        /// Vertex coefficient of mass matrix creation function to handle creation of new hexahedra
        void applyHexahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Hexahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Vertex coefficient of mass matrix destruction function to handle creation of new hexahedra
        void applyHexahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add hexahedron elements.
        virtual void ApplyTopologyChange(const core::topology::HexahedraAdded* /*event*/);
         /// Callback to remove hexahedron elements.
        virtual void ApplyTopologyChange(const core::topology::HexahedraRemoved* /*event*/);

    protected:
        MuscleActivation<DataTypes,TMassType>* m;
    };
    VertexMassHandler* vertexMassHandler;

    class EdgeMassHandler : public topology::TopologyDataHandler<Edge,MassVector>
    {
    public:
        EdgeMassHandler(MuscleActivation<DataTypes,TMassType>* _m, EdgeData<helper::vector<TMassType> >* _data) : topology::TopologyDataHandler<Edge,helper::vector<TMassType> >(_data), m(_m) {}

        /// Edge mass coefficient matrix creation function
        void applyCreateFunction(unsigned int edgeIndex, MassType & EdgeMass,
                const Edge&,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double >&);


        ///////////////////////// Functions on Triangles //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new triangles
        void applyTriangleCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Triangle >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new triangles
        void applyTriangleDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesAdded* /*event*/);
        /// Callback to remove triangles elements.
        void ApplyTopologyChange(const core::topology::TrianglesRemoved* /*event*/);


        ///////////////////////// Functions on Quads //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new quads
        void applyQuadCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Quad >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new quads
        void applyQuadDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add quads elements.
        void ApplyTopologyChange(const core::topology::QuadsAdded* /*event*/);
        /// Callback to remove quads elements.
        void ApplyTopologyChange(const core::topology::QuadsRemoved* /*event*/);


        ///////////////////////// Functions on Tetrahedron //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new tetrahedra
        void applyTetrahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Tetrahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new tetrahedra
        void applyTetrahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraAdded* /*event*/);
        /// Callback to remove tetrahedron elements.
        void ApplyTopologyChange(const core::topology::TetrahedraRemoved* /*event*/);


        ///////////////////////// Functions on Hexahedron //////////////////////////////////////

        /// Edge coefficient of mass matrix creation function to handle creation of new hexahedra
        void applyHexahedronCreation(const sofa::helper::vector< unsigned int >& /*indices*/,
                const sofa::helper::vector< Hexahedron >& /*elems*/,
                const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
                const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/);

        /// Edge coefficient of mass matrix destruction function to handle creation of new hexahedra
        void applyHexahedronDestruction(const sofa::helper::vector<unsigned int> & /*indices*/);

        /// Callback to add hexahedron elements.
        void ApplyTopologyChange(const core::topology::HexahedraAdded* /*event*/);
         /// Callback to remove hexahedron elements.
        void ApplyTopologyChange(const core::topology::HexahedraRemoved* /*event*/);

    protected:
        MuscleActivation<DataTypes,TMassType>* m;

    };

    EdgeMassHandler* edgeMassHandler;

    class TetrahedronMassHandler : public topology::TopologyDataHandler<Tetrahedron,MassVectorVector>
    {
    public:
        typedef typename DataTypes::Real Real;
        TetrahedronMassHandler(MuscleActivation<DataTypes,TMassType>* _m, TetrahedronData<helper::vector<MassVector> >* _data) : topology::TopologyDataHandler<Tetrahedron,helper::vector<MassVector> >(_data), m(_m) {}

        /// Edge mass coefficient matrix creation function
        void applyCreateFunction(unsigned int tetrahedronIndex, MassVector & tetrahedronMass,
                const Tetrahedron&,
                const sofa::helper::vector< unsigned int > &,
                const sofa::helper::vector< double >&);

               /// Edge coefficient of mass matrix destruction function to handle creation of new tetrahedra
        void applyDestructionFunction(const sofa::helper::vector<unsigned int> & /*indices*/);

    protected:
        MuscleActivation<DataTypes,TMassType>* m;
    };

    TetrahedronMassHandler* tetrahedronMassHandler;

};

} // namespace mass

} // namespace component

} // namespace sofa

#endif

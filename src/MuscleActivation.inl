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
#ifndef SOFA_COMPONENT_MASS_MuscleActivation_INL
#define SOFA_COMPONENT_MASS_MuscleActivation_INL

#include "MuscleActivation.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/template.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <SofaBaseTopology/TopologyData.inl>
#include <SofaBaseTopology/RegularGridTopology.h>
#include <SofaBaseMechanics/AddMToMatrixFunctor.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/vector.h>
#include <SofaBaseTopology/CommonAlgorithms.h>
#include <SofaBaseTopology/BezierTetrahedronSetGeometryAlgorithms.h>
#include <SofaBaseTopology/EdgeSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TetrahedronSetGeometryAlgorithms.h>
#include <SofaBaseTopology/QuadSetGeometryAlgorithms.h>
#include <SofaBaseTopology/HexahedronSetGeometryAlgorithms.h>
#include <SofaBaseMechanics/MechanicalObject.h>

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/rapidjson.h"
#include <iostream>
using namespace rapidjson;

#ifdef SOFA_SUPPORT_MOVING_FRAMES
#include <sofa/core/behavior/InertiaForce.h>
#endif

#ifdef SOFA_HAVE_EIGEN2
#include <Eigen/src/Core/util/ForwardDeclarations.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/LU/Inverse.h>
#include <Eigen/Geometry>
#endif

using namespace std;

namespace sofa
{

namespace component
{

namespace mass
{

template <class DataTypes, class MassType>
MuscleActivation<DataTypes, MassType>::MuscleActivation()
    : vertexMassInfo( initData(&vertexMassInfo, "vertexMass", "values of the particles masses on vertices") )
    , edgeMassInfo( initData(&edgeMassInfo, "edgeMass", "values of the particles masses on edges") )
    , tetrahedronMassInfo( initData(&tetrahedronMassInfo, "tetrahedronMass", "values of the particles masses for all control points inside a Bezier tetrahedron") )
    , m_massDensity( initData(&m_massDensity, (Real)1.0,"massDensity", "mass density that allows to compute the  particles masses from a mesh topology and geometry.\nOnly used if > 0") )
    , showCenterOfGravity( initData(&showCenterOfGravity, false, "showGravityCenter", "display the center of gravity of the system" ) )
    , showAxisSize( initData(&showAxisSize, (Real)1.0, "showAxisSizeFactor", "factor length of the axis displayed (only used for rigids)" ) )
    , lumping( initData(&lumping, true, "lumping","boolean if you need to use a lumped mass matrix") )
    , printMass( initData(&printMass, false, "printMass","boolean if you want to get the totalMass") )
    , f_graph( initData(&f_graph,"graph","Graph of the controlled potential") )
    , numericalIntegrationOrder( initData(&numericalIntegrationOrder,(size_t)2,"integrationOrder","The order of integration for numerical integration"))
    , numericalIntegrationMethod( initData(&numericalIntegrationMethod,(size_t)0,"numericalIntegrationMethod","The type of numerical integration method chosen"))
    , d_integrationMethod( initData(&d_integrationMethod,std::string("analytical"),"integrationMethod","\"exact\" if closed form expression for high order elements, \"analytical\" if closed form expression for affine element, \"numerical\" if numerical integration is chosen"))
	, f_fileNodes(initData(&f_fileNodes,"f_fileNodes","file with node indices that receive forces"))
	, f_fileActivation(initData(&f_fileActivation,"activation_order_file","order of nodes for muscle activation from IM"))
	, f_fileExtremes(initData(&f_fileExtremes,"file_extremes","file with node indices of muscle extremes"))
	, f_signalFile(initData(&f_signalFile,"signal_file","signals file"))
	, f_musclePositionFile(initData(&f_musclePositionFile,"muscle_position_file","file that matches muscles in alphabetical order with positiones muscles"))
	, f_activationConstant1(initData(&f_activationConstant1,(Real)10.0,"activation_scale1","scale the activation forces"))
	, f_activationConstant2(initData(&f_activationConstant2,(Real)10.0,"activation_scale2","scale the activation forces"))
	, f_activationPeriod1(initData(&f_activationPeriod1,(Real)1.0,"activation_period1","period of the activation forces 1"))
	, f_activationPeriod2(initData(&f_activationPeriod2,(Real)1.0,"activation_period2","period of the activation forces 2"))
	, f_maxActivatedMuscles(initData(&f_maxActivatedMuscles, 12, "max_activated_muscles","from 0 to 12, how many muscles will be activated"))
	, f_minActivatedMuscles(initData(&f_minActivatedMuscles, 0, "min_activated_muscles","from 0 to 12, how many muscles will be activated"))
	, f_neuronTimestep(initData(&f_neuronTimestep,(Real)0.001,"neuron_timestep","timestep of the neuron network in the file"))
	, f_activationStop1(initData(&f_activationStop1,(Real)0.001,"activation_stop1","activation_stop1"))
	, f_activationStop2(initData(&f_activationStop2,(Real)0.001,"activation_stop2","activation_stop2"))
	, f_correctiveForces( initData(&f_correctiveForces, false, "f_correctiveForces","f_correctiveForces") )
{
    f_graph.setWidget("graph");
	f_fileNodes.setGroup("Forces");
	f_fileActivation.setGroup("Forces");
	f_musclePositionFile.setGroup("Forces");
	f_activationConstant1.setGroup("Forces");
	f_activationConstant2.setGroup("Forces");
	f_activationPeriod1.setGroup("Forces");
	f_activationPeriod2.setGroup("Forces");
	f_minActivatedMuscles.setGroup("Forces");
	f_maxActivatedMuscles.setGroup("Forces");
	f_signalFile.setGroup("Forces");
	f_fileExtremes.setGroup("Forces");
	f_neuronTimestep.setGroup("Forces");
	f_correctiveForces.setGroup("Forces");

}

template <class DataTypes, class MassType>
MuscleActivation<DataTypes, MassType>::~MuscleActivation()
{
    if (vertexMassHandler) delete vertexMassHandler;
    if (edgeMassHandler) delete edgeMassHandler;
	if (tetrahedronMassHandler) delete tetrahedronMassHandler;
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyCreateFunction(unsigned int, MassType & VertexMass,
        const sofa::helper::vector< unsigned int > &,
        const sofa::helper::vector< double >&)
{
    VertexMass = 0;
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyCreateFunction(unsigned int, MassType & EdgeMass,
        const topology::Edge&,
        const sofa::helper::vector< unsigned int > &,
        const sofa::helper::vector< double >&)
{
    EdgeMass = 0;
}
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::TetrahedronMassHandler::applyCreateFunction(unsigned int tetra, MassVector & TetrahedronMass,
        const topology::Tetrahedron&,
        const sofa::helper::vector< unsigned int > &,
        const sofa::helper::vector< double >&)
{
	MuscleActivation<DataTypes, MassType> *MMM = this->m;

	if (MMM && (MMM->bezierTetraGeo) && (MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_BEZIERTETRAHEDRONSET))
	{	
		Real densityM = MMM->getMassDensity();
		topology::BezierDegreeType degree=MMM->bezierTetraGeo->getTopologyContainer()->getDegree();
		size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
		size_t nbMassEntries=nbControlPoints*(nbControlPoints+1)/2;

		if (TetrahedronMass.size()!=nbMassEntries) {
			TetrahedronMass.resize(nbMassEntries);
		}
		// set array to zero
		std::fill(TetrahedronMass.begin(),TetrahedronMass.end(),(MassType)0);
		sofa::helper::vector<MassType> lumpedVertexMass;
		lumpedVertexMass.resize(nbControlPoints);
		size_t i,j,k,rank;
		topology::VecPointID indexArray;
		/// get the global index of each control point in the tetrahedron
		MMM->bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfBezierPointsInTetrahedron(tetra,indexArray);
		std::fill(lumpedVertexMass.begin(),lumpedVertexMass.end(),(MassType)0);
		if (MMM->integrationMethod==MuscleActivation<DataTypes, MassType>::NUMERICAL_INTEGRATION) {
			sofa::helper::vector<Real> shapeFunctionValue;
			shapeFunctionValue.resize(nbControlPoints);
			// set array to zero

			/// get value of integration points
			topology::NumericalIntegrationDescriptor<Real,4> &nid=MMM->bezierTetraGeo->getTetrahedronNumericalIntegrationDescriptor();
            typename topology::NumericalIntegrationDescriptor<Real,4>::QuadraturePointArray qpa=nid.getQuadratureMethod((typename topology::NumericalIntegrationDescriptor<Real,4>::QuadratureMethod)MMM->numericalIntegrationMethod.getValue(),
				MMM->numericalIntegrationOrder.getValue());

			sofa::defaulttype::Vec<4,Real> bc;
			sofa::helper::vector<topology::TetrahedronBezierIndex> tbi=MMM->bezierTetraGeo->getTopologyContainer()->getTetrahedronBezierIndexArray();
			typename DataTypes::Real jac,weight;
			MassType tmpMass;

			// loop through the integration points
			for (i=0;i<qpa.size();++i) {
                typename topology::NumericalIntegrationDescriptor<Real,4>::QuadraturePoint qp=qpa[i];
				// the barycentric coordinate
				bc=qp.first;
				// the weight of the integration point
				weight=qp.second;
				// the Jacobian Derterminant of the integration point 
				jac=MMM->bezierTetraGeo->computeJacobian(tetra,bc)*densityM;
				/// prestore the shape function value for that integration point.
				for (j=0;j<nbControlPoints;j++) {
					shapeFunctionValue[j]=MMM->bezierTetraGeo->computeBernsteinPolynomial(tbi[j],bc);
				}
				// now loop through each pair of control point to compute the mass
				rank=0;
				for (j=0;j<nbControlPoints;j++) {
					// use the fact that the shapefunctions sum to 1 to get the lumped value
					lumpedVertexMass[j]+=shapeFunctionValue[j]*fabs(jac)*weight;

					for (k=j;k<nbControlPoints;k++,rank++) {
						/// compute the mass as the integral of product of the 2 shape functions multiplied by the Jacobian
						tmpMass=shapeFunctionValue[k]*shapeFunctionValue[j]*fabs(jac)*weight;
						TetrahedronMass[rank]+=tmpMass;
					}
				}
			}
	//		std::cerr<<"Mass Matrix= "<<TetrahedronMass<<std::endl;
	//		std::cerr<<"Lumped Mass Matrix= "<<lumpedVertexMass<<std::endl;
			// now updates the the mass matrix on each vertex.
			vector<MassType>& my_vertexMassInfo = *MMM->vertexMassInfo.beginEdit();
			for (j=0;j<nbControlPoints;j++) {
				my_vertexMassInfo[indexArray[j]]+=lumpedVertexMass[j];
			}
		} else if ((MMM->integrationMethod==MuscleActivation<DataTypes, MassType>::AFFINE_ELEMENT_INTEGRATION) || 
            (MMM->bezierTetraGeo->isBezierTetrahedronAffine(tetra,(MMM->bezierTetraGeo->getDOF()->read(core::ConstVecCoordId::restPosition())->getValue()) )))
		{
			/// affine mass simple computation
			
			Real totalMass= densityM*MMM->tetraGeo->computeRestTetrahedronVolume(tetra);
			Real mass=totalMass/(topology::binomial<typename DataTypes::Real>(degree,degree)*topology::binomial<typename DataTypes::Real>(2*degree,3));
			sofa::helper::vector<topology::TetrahedronBezierIndex> tbiArray;
			topology::TetrahedronBezierIndex tbi1,tbi2;
			tbiArray=MMM->bezierTetraGeo->getTopologyContainer()->getTetrahedronBezierIndexArray();
			rank=0;
			for (j=0;j<nbControlPoints;j++) {
				tbi1=tbiArray[j];
				for (k=j;k<nbControlPoints;k++,rank++) {
					tbi2=tbiArray[k];
					TetrahedronMass[rank]+=mass*topology::binomialVector<4,typename DataTypes::Real>(tbi1,tbi2);
	//				std::cerr<<" tbi = "<<tbi1<<" "<<tbi2<<" ="<<TetrahedronMass[rank]<<std::endl;
				}
			}
			// mass for mass lumping
			mass=totalMass/nbControlPoints;
			// now updates the the mass matrix on each vertex.
			vector<MassType>& my_vertexMassInfo = *MMM->vertexMassInfo.beginEdit();
			for (j=0;j<nbControlPoints;j++) {
				my_vertexMassInfo[indexArray[j]]+=mass;
			}
		} else {
			/// exact computation
			sofa::helper::vector<topology::TetrahedronBezierIndex> tbiArray,tbiDerivArray,multinomialArray;
			sofa::helper::vector<unsigned char> multinomialScalarArray;
			/// use the rest configuration
//			const typename DataTypes::VecCoord &p=(MMM->bezierTetraGeo->getDOF()->read(core::ConstVecCoordId::restPosition())->getValue());

			Real factor;
			MassType tmpMass;

			tbiArray=MMM->bezierTetraGeo->getTopologyContainer()->getTetrahedronBezierIndexArray();
			tbiDerivArray=MMM->bezierTetraGeo->getTopologyContainer()->getTetrahedronBezierIndexArrayOfGivenDegree(degree-1);
			sofa::helper::vector<topology::LocalTetrahedronIndex> correspondanceArray=
				MMM->bezierTetraGeo->getTopologyContainer()->getMapOfTetrahedronBezierIndexArrayFromInferiorDegree();
			typename DataTypes::Coord dp1,dp2,dp3,dpos,tmp;
			multinomialArray.resize(5);
			multinomialScalarArray.resize(5);
			multinomialScalarArray[0]=degree-1;
			multinomialScalarArray[1]=degree-1;
			multinomialScalarArray[2]=degree-1;
			multinomialScalarArray[3]=degree;
			multinomialScalarArray[4]=degree;
			size_t l,m;
			factor=6*topology::multinomial<Real>(5*degree-3,multinomialScalarArray)*topology::binomial<Real>(5*degree-3,3)/(degree*degree*degree*densityM);
			for (i=0;i<tbiDerivArray.size();++i) 
			{
				dp1=MMM->bezierTetraGeo->getPointRestPosition(indexArray[correspondanceArray[i][0]])-MMM->bezierTetraGeo->getPointRestPosition(indexArray[correspondanceArray[i][3]]);
				multinomialArray[0]=tbiDerivArray[i];
				for (j=0;j<tbiDerivArray.size();++j) 
				{
					dp2=MMM->bezierTetraGeo->getPointRestPosition(indexArray[correspondanceArray[j][1]])-
						MMM->bezierTetraGeo->getPointRestPosition(indexArray[correspondanceArray[j][3]]);
                    using topology::cross;
					tmp=cross<Real>(dp1,dp2);
					multinomialArray[1]=tbiDerivArray[j];
					rank=0;
					for (l=0;l<nbControlPoints;l++) {
						multinomialArray[3]=tbiArray[l];
						for (m=l;m<nbControlPoints;m++,rank++) {
							multinomialArray[4]=tbiArray[m];
							// set dp3 to 0 in a generic way
							std::fill(dp3.begin(),dp3.end(),(Real)0);
							for (k=0;k<tbiDerivArray.size();++k) 
							{
								multinomialArray[2]=tbiDerivArray[k];
								dpos=MMM->bezierTetraGeo->getPointRestPosition(indexArray[correspondanceArray[k][2]])-
									MMM->bezierTetraGeo->getPointRestPosition(indexArray[correspondanceArray[k][3]]);
								dpos*=topology::multinomialVector<4,Real>(multinomialArray);
								dp3+=dpos;
							}
							dp3/=factor;
							tmpMass=fabs(dp3*tmp);
							TetrahedronMass[rank]+=tmpMass;
							lumpedVertexMass[l]+=tmpMass;
							if (m>l)
								lumpedVertexMass[m]+=tmpMass;
						}

					}
				}

			}

			// now updates the the mass matrix on each vertex.
			vector<MassType>& my_vertexMassInfo = *MMM->vertexMassInfo.beginEdit();
			for (j=0;j<nbControlPoints;j++) {
				my_vertexMassInfo[indexArray[j]]+=lumpedVertexMass[j];
			}
		}
	}
}

// -------------------------------------------------------
// ------- Triangle Creation/Destruction functions -------
// -------------------------------------------------------
//{

/// Creation fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyTriangleCreation(const sofa::helper::vector< unsigned int >& triangleAdded,
        const sofa::helper::vector< topology::Triangle >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TRIANGLESET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<triangleAdded.size(); ++i)
        {
            // Get the triangle to be added
            const topology::Triangle &t = MMM->_topology->getTriangle(triangleAdded[i]);

            // Compute rest mass of conserne triangle = density * triangle surface.
            if(MMM->triangleGeo)
            {
                mass=(densityM * MMM->triangleGeo->computeRestTriangleArea(triangleAdded[i]))/(typename DataTypes::Real)6.0;
            }

            // Adding mass
            for (unsigned int j=0; j<3; ++j)
                VertexMasses[ t[j] ] += mass;
        }
    }
}

/// Creation fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyTriangleCreation(const sofa::helper::vector< unsigned int >& triangleAdded,
        const sofa::helper::vector< topology::Triangle >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TRIANGLESET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<triangleAdded.size(); ++i)
        {
            // Get the edgesInTriangle to be added
            const topology::EdgesInTriangle &te = MMM->_topology->getEdgesInTriangle(triangleAdded[i]);

            // Compute rest mass of conserne triangle = density * triangle surface.
            if(MMM->triangleGeo)
            {
                mass=(densityM * MMM->triangleGeo->computeRestTriangleArea(triangleAdded[i]))/(typename DataTypes::Real)12.0;
            }

            // Adding mass edges of concerne triangle
            for (unsigned int j=0; j<3; ++j)
                EdgeMasses[ te[j] ] += mass;
        }
    }
}


/// Destruction fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyTriangleDestruction(const sofa::helper::vector< unsigned int >& triangleRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TRIANGLESET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<triangleRemoved.size(); ++i)
        {
            // Get the triangle to be removed
            const topology::Triangle &t = MMM->_topology->getTriangle(triangleRemoved[i]);

            // Compute rest mass of conserne triangle = density * triangle surface.
            if(MMM->triangleGeo)
            {
                mass=(densityM * MMM->triangleGeo->computeRestTriangleArea(triangleRemoved[i]))/(typename DataTypes::Real)6.0;
            }

            // Removing mass
            for (unsigned int j=0; j<3; ++j)
                VertexMasses[ t[j] ] -= mass;
        }
    }
}


/// Destruction fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyTriangleDestruction(const sofa::helper::vector< unsigned int >& triangleRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TRIANGLESET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<triangleRemoved.size(); ++i)
        {
            // Get the triangle to be removed
            const topology::EdgesInTriangle &te = MMM->_topology->getEdgesInTriangle(triangleRemoved[i]);

            // Compute rest mass of conserne triangle = density * triangle surface.
            if(MMM->triangleGeo)
            {
                mass=(densityM * MMM->triangleGeo->computeRestTriangleArea(triangleRemoved[i]))/(typename DataTypes::Real)12.0;
            }

            // Removing mass edges of concerne triangle
            for (unsigned int j=0; j<3; ++j)
                EdgeMasses[ te[j] ] -= mass;
        }
    }
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::TrianglesAdded* e)
{
    const sofa::helper::vector<unsigned int> &triangleAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Triangle> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyTriangleCreation(triangleAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::TrianglesRemoved* e)
{
    const sofa::helper::vector<unsigned int> &triangleRemoved = e->getArray();

    applyTriangleDestruction(triangleRemoved);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::TrianglesAdded* e)
{
    const sofa::helper::vector<unsigned int> &triangleAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Triangle> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyTriangleCreation(triangleAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::TrianglesRemoved* e)
{
    const sofa::helper::vector<unsigned int> &triangleRemoved = e->getArray();

    applyTriangleDestruction(triangleRemoved);
}

// }

// ---------------------------------------------------
// ------- Quad Creation/Destruction functions -------
// ---------------------------------------------------
//{

/// Creation fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyQuadCreation(const sofa::helper::vector< unsigned int >& quadAdded,
        const sofa::helper::vector< topology::Quad >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_QUADSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<quadAdded.size(); ++i)
        {
            // Get the quad to be added
            const topology::Quad &q = MMM->_topology->getQuad(quadAdded[i]);

            // Compute rest mass of conserne quad = density * quad surface.
            if(MMM->quadGeo)
            {
                mass=(densityM * MMM->quadGeo->computeRestQuadArea(quadAdded[i]))/(typename DataTypes::Real)8.0;
            }

            // Adding mass
            for (unsigned int j=0; j<4; ++j)
                VertexMasses[ q[j] ] += mass;
        }
    }
}


/// Creation fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyQuadCreation(const sofa::helper::vector< unsigned int >& quadAdded,
        const sofa::helper::vector< topology::Quad >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_QUADSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<quadAdded.size(); ++i)
        {
            // Get the EdgesInQuad to be added
            const topology::EdgesInQuad &qe = MMM->_topology->getEdgesInQuad(quadAdded[i]);

            // Compute rest mass of conserne quad = density * quad surface.
            if(MMM->quadGeo)
            {
                mass=(densityM * MMM->quadGeo->computeRestQuadArea(quadAdded[i]))/(typename DataTypes::Real)16.0;
            }

            // Adding mass edges of concerne quad
            for (unsigned int j=0; j<4; ++j)
                EdgeMasses[ qe[j] ] += mass;
        }
    }
}


/// Destruction fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyQuadDestruction(const sofa::helper::vector< unsigned int >& quadRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_QUADSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<quadRemoved.size(); ++i)
        {
            // Get the quad to be removed
            const topology::Quad &q = MMM->_topology->getQuad(quadRemoved[i]);

            // Compute rest mass of conserne quad = density * quad surface.
            if(MMM->quadGeo)
            {
                mass=(densityM * MMM->quadGeo->computeRestQuadArea(quadRemoved[i]))/(typename DataTypes::Real)8.0;
            }

            // Removing mass
            for (unsigned int j=0; j<4; ++j)
                VertexMasses[ q[j] ] -= mass;
        }
    }
}

/// Destruction fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyQuadDestruction(const sofa::helper::vector< unsigned int >& quadRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_QUADSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<quadRemoved.size(); ++i)
        {
            // Get the EdgesInQuad to be removed
            const topology::EdgesInQuad &qe = MMM->_topology->getEdgesInQuad(quadRemoved[i]);

            // Compute rest mass of conserne quad = density * quad surface.
            if(MMM->quadGeo)
            {
                mass=(densityM * MMM->quadGeo->computeRestQuadArea(quadRemoved[i]))/(typename DataTypes::Real)16.0;
            }

            // Removing mass edges of concerne quad
            for (unsigned int j=0; j<4; ++j)
                EdgeMasses[ qe[j] ] -= mass/2;
        }
    }
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::QuadsAdded* e)
{
    const sofa::helper::vector<unsigned int> &quadAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Quad> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyQuadCreation(quadAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::QuadsRemoved* e)
{
    const sofa::helper::vector<unsigned int> &quadRemoved = e->getArray();

    applyQuadDestruction(quadRemoved);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::QuadsAdded* e)
{
    const sofa::helper::vector<unsigned int> &quadAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Quad> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyQuadCreation(quadAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::QuadsRemoved* e)
{
    const sofa::helper::vector<unsigned int> &quadRemoved = e->getArray();

    applyQuadDestruction(quadRemoved);
}

// }



// ----------------------------------------------------------
// ------- Tetrahedron Creation/Destruction functions -------
// ----------------------------------------------------------
//{

/// Creation fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyTetrahedronCreation(const sofa::helper::vector< unsigned int >& tetrahedronAdded,
        const sofa::helper::vector< topology::Tetrahedron >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TETRAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<tetrahedronAdded.size(); ++i)
        {
            // Get the tetrahedron to be added
            const topology::Tetrahedron &t = MMM->_topology->getTetrahedron(tetrahedronAdded[i]);

            // Compute rest mass of conserne tetrahedron = density * tetrahedron volume.
            if(MMM->tetraGeo)
            {
                mass=(densityM * MMM->tetraGeo->computeRestTetrahedronVolume(tetrahedronAdded[i]))/(typename DataTypes::Real)10.0;
            }

            // Adding mass
            for (unsigned int j=0; j<4; ++j)
                VertexMasses[ t[j] ] += mass;
        }
    }
}


/// Creation fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyTetrahedronCreation(const sofa::helper::vector< unsigned int >& tetrahedronAdded,
        const sofa::helper::vector< topology::Tetrahedron >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TETRAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<tetrahedronAdded.size(); ++i)
        {
            // Get the edgesInTetrahedron to be added
            const topology::EdgesInTetrahedron &te = MMM->_topology->getEdgesInTetrahedron(tetrahedronAdded[i]);

            // Compute rest mass of conserne triangle = density * tetrahedron volume.
            if(MMM->tetraGeo)
            {
                mass=(densityM * MMM->tetraGeo->computeRestTetrahedronVolume(tetrahedronAdded[i]))/(typename DataTypes::Real)20.0;
            }

            // Adding mass edges of concerne triangle
            for (unsigned int j=0; j<6; ++j)
                EdgeMasses[ te[j] ] += mass;
        }
    }
}


/// Destruction fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyTetrahedronDestruction(const sofa::helper::vector< unsigned int >& tetrahedronRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TETRAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<tetrahedronRemoved.size(); ++i)
        {
            // Get the tetrahedron to be removed
            const topology::Tetrahedron &t = MMM->_topology->getTetrahedron(tetrahedronRemoved[i]);

            // Compute rest mass of conserne tetrahedron = density * tetrahedron volume.
            if(MMM->tetraGeo)
            {
                mass=(densityM * MMM->tetraGeo->computeRestTetrahedronVolume(tetrahedronRemoved[i]))/(typename DataTypes::Real)10.0;
            }

            // Removing mass
            for (unsigned int j=0; j<4; ++j)
                VertexMasses[ t[j] ] -= mass;
        }
    }
}

/// Destruction fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyTetrahedronDestruction(const sofa::helper::vector< unsigned int >& tetrahedronRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_TETRAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<tetrahedronRemoved.size(); ++i)
        {
            // Get the edgesInTetrahedron to be removed
            const topology::EdgesInTetrahedron &te = MMM->_topology->getEdgesInTetrahedron(tetrahedronRemoved[i]);

            // Compute rest mass of conserne triangle = density * tetrahedron volume.
            if(MMM->tetraGeo)
            {
                mass=(densityM * MMM->tetraGeo->computeRestTetrahedronVolume(tetrahedronRemoved[i]))/(typename DataTypes::Real)20.0;
            }

            // Removing mass edges of concerne triangle
            for (unsigned int j=0; j<6; ++j)
                EdgeMasses[ te[j] ] -= mass; //?
        }
    }
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::TetrahedraAdded* e)
{
    const sofa::helper::vector<unsigned int> &tetraAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Tetrahedron> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyTetrahedronCreation(tetraAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::TetrahedraRemoved* e)
{
    const sofa::helper::vector<unsigned int> &tetraRemoved = e->getArray();

    applyTetrahedronDestruction(tetraRemoved);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::TetrahedraAdded* e)
{
    const sofa::helper::vector<unsigned int> &tetraAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Tetrahedron> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyTetrahedronCreation(tetraAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::TetrahedraRemoved* e)
{
    const sofa::helper::vector<unsigned int> &tetraRemoved = e->getArray();

    applyTetrahedronDestruction(tetraRemoved);
}

// }


// ---------------------------------------------------------
// ------- Hexahedron Creation/Destruction functions -------
// ---------------------------------------------------------
//{

/// Creation fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyHexahedronCreation(const sofa::helper::vector< unsigned int >& hexahedronAdded,
        const sofa::helper::vector< topology::Hexahedron >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_HEXAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<hexahedronAdded.size(); ++i)
        {
            // Get the hexahedron to be added
            const topology::Hexahedron &h = MMM->_topology->getHexahedron(hexahedronAdded[i]);

            // Compute rest mass of conserne hexahedron = density * hexahedron volume.
            if(MMM->hexaGeo)
            {
                mass=(densityM * MMM->hexaGeo->computeRestHexahedronVolume(hexahedronAdded[i]))/(typename DataTypes::Real)20.0;
            }

            // Adding mass
            for (unsigned int j=0; j<8; ++j)
                VertexMasses[ h[j] ] += mass;
        }
    }
}


/// Creation fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyHexahedronCreation(const sofa::helper::vector< unsigned int >& hexahedronAdded,
        const sofa::helper::vector< topology::Hexahedron >& /*elems*/,
        const sofa::helper::vector< sofa::helper::vector< unsigned int > >& /*ancestors*/,
        const sofa::helper::vector< sofa::helper::vector< double > >& /*coefs*/)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;

    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_HEXAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<hexahedronAdded.size(); ++i)
        {
            // Get the EdgesInHexahedron to be added
            const topology::EdgesInHexahedron &he = MMM->_topology->getEdgesInHexahedron(hexahedronAdded[i]);

            // Compute rest mass of conserne hexahedron = density * hexahedron volume.
            if(MMM->hexaGeo)
            {
                mass=(densityM * MMM->hexaGeo->computeRestHexahedronVolume(hexahedronAdded[i]))/(typename DataTypes::Real)40.0;
            }

            // Adding mass edges of concerne triangle
            for (unsigned int j=0; j<12; ++j)
                EdgeMasses[ he[j] ] += mass;
        }
    }
}


/// Destruction fonction for mass stored on vertices
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::applyHexahedronDestruction(const sofa::helper::vector< unsigned int >& hexahedronRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_HEXAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > VertexMasses ( MMM->vertexMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<hexahedronRemoved.size(); ++i)
        {
            // Get the hexahedron to be removed
            const topology::Hexahedron &h = MMM->_topology->getHexahedron(hexahedronRemoved[i]);

            // Compute rest mass of conserne hexahedron = density * hexahedron volume.
            if(MMM->hexaGeo)
            {
                mass=(densityM * MMM->hexaGeo->computeRestHexahedronVolume(hexahedronRemoved[i]))/(typename DataTypes::Real)20.0;
            }

            // Removing mass
            for (unsigned int j=0; j<8; ++j)
                VertexMasses[ h[j] ] -= mass;
        }
    }
}


/// Destruction fonction for mass stored on edges
template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::applyHexahedronDestruction(const sofa::helper::vector< unsigned int >& hexahedronRemoved)
{
    MuscleActivation<DataTypes, MassType> *MMM = this->m;
    if (MMM && MMM->getMassTopologyType()==MuscleActivation<DataTypes, MassType>::TOPOLOGY_HEXAHEDRONSET)
    {
        helper::WriteAccessor< Data< vector<MassType> > > EdgeMasses ( MMM->edgeMassInfo );
        // Initialisation
        typename DataTypes::Real densityM = MMM->getMassDensity();
        typename DataTypes::Real mass = (typename DataTypes::Real) 0;

        for (unsigned int i = 0; i<hexahedronRemoved.size(); ++i)
        {
            // Get the EdgesInHexahedron to be removed
            const topology::EdgesInHexahedron &he = MMM->_topology->getEdgesInHexahedron(hexahedronRemoved[i]);

            // Compute rest mass of conserne hexahedron = density * hexahedron volume.
            if(MMM->hexaGeo)
            {
                mass=(densityM * MMM->hexaGeo->computeRestHexahedronVolume(hexahedronRemoved[i]))/(typename DataTypes::Real)40.0;
            }

            // Removing mass edges of concerne triangle
            for (unsigned int j=0; j<12; ++j)
                EdgeMasses[ he[j] ] -= mass;
        }
    }
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::HexahedraAdded* e)
{
    const sofa::helper::vector<unsigned int> &hexaAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Hexahedron> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyHexahedronCreation(hexaAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::VertexMassHandler::ApplyTopologyChange(const core::topology::HexahedraRemoved* e)
{
    const sofa::helper::vector<unsigned int> &hexaRemoved = e->getArray();

    applyHexahedronDestruction(hexaRemoved);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::HexahedraAdded* e)
{
    const sofa::helper::vector<unsigned int> &hexaAdded = e->getIndexArray();
    const sofa::helper::vector<topology::Hexahedron> &elems = e->getElementArray();
    const sofa::helper::vector<sofa::helper::vector<unsigned int> > & ancestors = e->ancestorsList;
    const sofa::helper::vector<sofa::helper::vector<double> > & coefs = e->coefs;

    applyHexahedronCreation(hexaAdded, elems, ancestors, coefs);
}

template< class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::EdgeMassHandler::ApplyTopologyChange(const core::topology::HexahedraRemoved* e)
{
    const sofa::helper::vector<unsigned int> &hexaRemoved = e->getArray();

    applyHexahedronDestruction(hexaRemoved);
}

// }



template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::init()
{
	std::cout << "init Mesh" << std::endl;

    this->Inherited::init();
    massLumpingCoeff = 0.0;

	if (d_integrationMethod.getValue() == "analytical")
		integrationMethod= AFFINE_ELEMENT_INTEGRATION;
	else if (d_integrationMethod.getValue() == "numerical") 
		integrationMethod= NUMERICAL_INTEGRATION;
	else if (d_integrationMethod.getValue() == "exact") 
		integrationMethod= EXACT_INTEGRATION;
	else
	{
		serr << "cannot recognize method "<< d_integrationMethod.getValue() << ". Must be either  \"exact\", \"analytical\" or \"numerical\"" << sendl;
	}

    _topology = this->getContext()->getMeshTopology();
    savedMass = m_massDensity.getValue();

    //    sofa::core::objectmodel::Tag mechanicalTag(m_tagMeshMechanics.getValue());
    //    this->getContext()->get(triangleGeo, mechanicalTag,sofa::core::objectmodel::BaseContext::SearchUp);

    this->getContext()->get(edgeGeo);
    this->getContext()->get(triangleGeo);
    this->getContext()->get(quadGeo);
    this->getContext()->get(tetraGeo);
    this->getContext()->get(hexaGeo);
	this->getContext()->get(bezierTetraGeo);

    // add the functions to handle topology changes for Vertex informations
    vertexMassHandler = new VertexMassHandler(this, &vertexMassInfo);
    vertexMassInfo.createTopologicalEngine(_topology, vertexMassHandler);
    vertexMassInfo.linkToEdgeDataArray();
    vertexMassInfo.linkToTriangleDataArray();
    vertexMassInfo.linkToQuadDataArray();
    vertexMassInfo.linkToTetrahedronDataArray();
    vertexMassInfo.linkToHexahedronDataArray();
    vertexMassInfo.registerTopologicalData();

    // add the functions to handle topology changes for Edge informations
    edgeMassHandler = new EdgeMassHandler(this, &edgeMassInfo);
    edgeMassInfo.createTopologicalEngine(_topology, edgeMassHandler);
    edgeMassInfo.linkToTriangleDataArray();
    edgeMassInfo.linkToQuadDataArray();
    edgeMassInfo.linkToTetrahedronDataArray();
    edgeMassInfo.linkToHexahedronDataArray();
    edgeMassInfo.registerTopologicalData();

	if (bezierTetraGeo) {
		// for Bezier Tetrahedra add the functions to handle topology changes for Tetrahedron informations
		tetrahedronMassHandler = new TetrahedronMassHandler(this, &tetrahedronMassInfo);
		tetrahedronMassInfo.createTopologicalEngine(_topology, tetrahedronMassHandler);
		tetrahedronMassInfo.linkToTetrahedronDataArray();

	}

    if ((vertexMassInfo.getValue().size()==0 || edgeMassInfo.getValue().size()==0) && (_topology!=0))
        reinit();

    //Reset the graph
    f_graph.beginEdit()->clear();
    f_graph.endEdit();

    this->copyVertexMass();

	//load the file that contains the indices of the nodes that will receive the forces
	if (!f_fileNodes.getValue().empty())
	{
		this->loadForceNodeFile(f_fileNodes.getFullPath().c_str());
	}

	//load the file that contains the indices of the nodes in order regarding muscles
	if (!f_fileActivation.getValue().empty())
	{
		this->loadActivationOrderFile(f_fileActivation.getFullPath().c_str());
	}

	std::cout << "loadMuscleExtremes" << std::endl;
	if (!f_fileExtremes.getValue().empty())
	{
		this->loadMuscleExtremes(f_fileExtremes.getFullPath().c_str());
	}

	std::cout << "load input signals" << std::endl;
	if (!f_signalFile.getValue().empty())
	{
		this->loadSignalFile(f_signalFile.getFullPath().c_str());
	}

	std::cout << "load muscle order" << std::endl;
	if (!f_musclePositionFile.getValue().empty())
	{
		this->loadMusclePositions(f_musclePositionFile.getFullPath().c_str());
	}

	m_activationConstant1 = f_activationConstant1.getValue();
	m_activationConstant2 = f_activationConstant2.getValue();
	m_activationPeriod1 = f_activationPeriod1.getValue();
	m_activationPeriod2 = f_activationPeriod2.getValue();
	m_activationStop1 = f_activationStop1.getValue();
	m_activationStop2 = f_activationStop2.getValue();
	m_minActivatedMuscles = f_minActivatedMuscles.getValue();
	m_maxActivatedMuscles = f_maxActivatedMuscles.getValue();

	if(m_minActivatedMuscles < 0) m_minActivatedMuscles = 0;

	for(unsigned int i = 0; i < 8; i++){
		std::vector<double> aux;
		for(unsigned int j = 0; j < 12; j++){
			aux.push_back(0.0);
		}
		m_saveSignals.push_back(aux);
	}

	int min = 0;
	//initialize forceHelper
	for(unsigned int i = 0; i < m_orderedMuscles.size(); i++){
		for(unsigned int j = m_minActivatedMuscles; j < m_orderedMuscles.at(i).size() && j < m_maxActivatedMuscles; j++){
			for(unsigned int k = 0; k < m_orderedMuscles.at(i).at(j).size(); k++){
				for(unsigned int l = 0; l < m_orderedMuscles.at(i).at(j).at(k).size(); l++){
					if(min < m_orderedMuscles[i][j][k][l]){
						min = m_orderedMuscles[i][j][k][l];
					}
				}
			}
		}
	}

	for(unsigned int i = 0; i < min; i++){
		Eigen::Vector3d force;
		force(0)=force(1)=force(2)= 0;
		m_forceHelper.push_back(force);
	}
}

template <class DataTypes, class MassType>
bool MuscleActivation<DataTypes, MassType>::loadMuscleExtremes(const char* filename){
	/*std::string meshFilename(filename);
	if (!sofa::helper::system::DataRepository.findFile (meshFilename))
	{
		serr << "Mesh \""<< filename <<"\" not found"<< sendl;
		return false;
	}
	this->f_fileExtremes.setValue( filename );*/
	std::ifstream file(filename, std::ios_base::in);
	if (!file.good()) return false;

	std::string line;			
	std::vector<unsigned int> aux;
	while(std::getline(file, line)){
		if(line == "#"){
			m_extremeIndices.push_back(aux);
			aux.clear();
			std::getline(file, line);
		}
		int aux2 = std::stoi(line);
		aux.push_back(aux2);
	}
	m_extremeIndices.push_back(aux);

	//check file read
	/*for(unsigned int i = 0; i < m_extremeIndices.size(); i++){
		for(unsigned int j = 0; j < m_extremeIndices.at(i).size(); j++){
			std::cout << i << "," << j << ": " << m_extremeIndices.at(i).at(j) << std::endl;
		}
	}*/

	for(unsigned int i = 0; i < m_extremeIndices.size(); i++){
		std::vector<std::vector<std::vector<unsigned int>>> aux3;
		for(unsigned int j = 0; j < m_extremeIndices.at(i).size() - 2; j = j+2){
			int ind = m_extremeIndices.at(i).at(j) - ((m_extremeIndices.at(i).at(j)-330)%8);
			int ind2 = m_extremeIndices.at(i).at(j+2) - ((m_extremeIndices.at(i).at(j+2)-330)%8);
			std::vector<std::vector<unsigned int>> aux2;
			while(ind <= ind2){

				//four corners of each row in the muscle
				double maxY = -999999999;
				double maxZ = -999999999;
				double minY = 999999999;
				double minZ = 999999999;
				unsigned int maxYIndex = 0;
				unsigned int maxZIndex = 0;
				unsigned int minYIndex = 0;
				unsigned int minZIndex = 0;
				unsigned int Y1Z1Index = 0;
				unsigned int Y1Z2Index = 0;
				unsigned int Y2Z1Index = 0;
				unsigned int Y2Z2Index = 0;
				for(unsigned int k = 0; k < 8; k++){
					if(_topology->getPY(ind+k) > maxY) {
						maxYIndex = ind+k;
						maxY = _topology->getPY(ind+k);
					}
					if(_topology->getPZ(ind+k) > maxZ){
						maxZIndex = ind+k;
						maxZ = _topology->getPZ(ind+k);
					}
					if(_topology->getPY(ind+k) < minY){
						minYIndex = ind+k;
						minY = _topology->getPY(ind+k);
					}
					if(_topology->getPZ(ind+k) < minZ){
						minZIndex = ind+k;
						minZ = _topology->getPZ(ind+k);
					}
				}

				//order other four nodes
				std::vector<unsigned int> remaining;
				for(unsigned int k = 0; k < 8; k++){
					if(ind + k == maxYIndex) continue;
					if(ind + k == maxZIndex) continue;
					if(ind + k == minYIndex) continue;
					if(ind + k == minZIndex) continue;
					bool A = true;
					for(unsigned int l = 0; l < remaining.size(); l++){
						if(_topology->getPY(ind+k) < _topology->getPY(remaining.at(l))){
							remaining.insert(remaining.begin()+l,ind+k);
							A =false;
							break;
						}
					}
					if(A) remaining.push_back(ind + k);
				}

				if(_topology->getPZ(remaining.at(0)) > _topology->getPZ(remaining.at(1))){
					Y1Z1Index = remaining.at(1);
					Y1Z2Index = remaining.at(0);
				}
				else{
					Y1Z1Index = remaining.at(0);
					Y1Z2Index = remaining.at(1);
				}

				if(_topology->getPZ(remaining.at(2)) > _topology->getPZ(remaining.at(3))){
					Y2Z1Index = remaining.at(3);
					Y2Z2Index = remaining.at(2);
				}
				else{
					Y2Z1Index = remaining.at(2);
					Y2Z2Index = remaining.at(3);
				}

				//check
				if(maxYIndex == maxZIndex) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxYIndex == minYIndex) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxYIndex == minZIndex) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxYIndex == Y1Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxYIndex == Y1Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxYIndex == Y2Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxYIndex == Y2Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxZIndex == minYIndex) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxZIndex == minZIndex) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxZIndex == Y1Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxZIndex == Y1Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxZIndex == Y2Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(maxZIndex == Y2Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minYIndex == minZIndex) std::cout << "ERROR: repeated indices" << std::endl;
				if(minYIndex == Y1Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minYIndex == Y1Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minYIndex == Y2Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minYIndex == Y2Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minZIndex == Y1Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minZIndex == Y1Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minZIndex == Y2Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(minZIndex == Y2Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(Y1Z1Index == Y1Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(Y1Z1Index == Y2Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(Y1Z1Index == Y2Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(Y1Z2Index == Y2Z1Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(Y1Z2Index == Y2Z2Index) std::cout << "ERROR: repeated indices" << std::endl;
				if(Y2Z1Index == Y2Z2Index) std::cout << "ERROR: repeated indices" << std::endl;

				//insert all the nodes in order
				std::vector<unsigned int> aux1;
				aux1.push_back(minYIndex);
				aux1.push_back(Y1Z2Index);
				aux1.push_back(maxZIndex);
				aux1.push_back(Y2Z2Index);
				aux1.push_back(maxYIndex);
				aux1.push_back(Y2Z1Index);
				aux1.push_back(minZIndex);
				aux1.push_back(Y1Z1Index);
				aux2.push_back(aux1);
				
				//next row in the muscle
				ind = ind + 8;
			}
			aux3.push_back(aux2);
		}
		m_orderedMuscles.push_back(aux3);
	}

	return true;
}

template <class DataTypes, class MassType>
std::vector<std::string> &MuscleActivation<DataTypes, MassType>::split(std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

template <class DataTypes, class MassType>
std::vector<std::string> MuscleActivation<DataTypes, MassType>::split(std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

template <class DataTypes, class MassType>
bool MuscleActivation<DataTypes, MassType>::loadActivationOrderFile(const char* filename){
	/*std::string meshFilename(filename);
	if (!sofa::helper::system::DataRepository.findFile (meshFilename))
	{
		serr << "Mesh \""<< filename <<"\" not found"<< sendl;
		return false;
	}
	this->f_fileNodes.setValue( filename );*/

	std::ifstream file(filename);
	if (!file.good()) return false;

	std::string line;
				

	while(std::getline(file, line)){
		std::vector<std::string> values = split(line,',');
		std::vector<unsigned int> intValues;
		for(unsigned int i = 0; i < values.size(); i++){
			intValues.push_back(std::stoi(values.at(i)));
		}
		m_activationOrder.push_back(intValues);
	}

	//check
	/*for(unsigned int i = 0; i < m_activationOrder.size(); i++){
		std::cout << i << ": ";
		for(unsigned int j = 0; j < m_activationOrder.at(i).size(); j++){
			std::cout << m_activationOrder.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}*/

	return true;
}

template <class DataTypes, class MassType>
bool MuscleActivation<DataTypes, MassType>::loadSignalFile(const char* filename){
	/*std::string meshFilename(filename);
	if (!sofa::helper::system::DataRepository.findFile (meshFilename))
	{
		serr << "Mesh \""<< filename <<"\" not found"<< sendl;
		return false;
	}
	this->f_fileNodes.setValue( filename );*/

	std::ifstream file(filename);
	if (!file.good()) return false;

	std::string line;
				

	while(std::getline(file, line)){
		std::vector<std::string> values = split(line,',');
		std::vector<double> doubleValues;
		for(unsigned int i = 0; i < values.size(); i++){
			doubleValues.push_back(std::stof(values.at(i)));
		}
		m_activationSignals.push_back(doubleValues);
	}

	//check
	/*for(unsigned int i = 0; i < m_activationSignals.size(); i++){
		std::cout << i << ": ";
		for(unsigned int j = 0; j < m_activationSignals.at(i).size(); j++){
			std::cout << m_activationSignals.at(i).at(j) << " ";
		}
		std::cout << std::endl;
	}
	*/
	return true;
}

template <class DataTypes, class MassType>
bool MuscleActivation<DataTypes, MassType>::loadMusclePositions(const char* filename){
	/*std::string meshFilename(filename);
	if (!sofa::helper::system::DataRepository.findFile (meshFilename))
	{
		serr << "Mesh \""<< filename <<"\" not found"<< sendl;
		return false;
	}
	this->f_fileMuscleExtremes.setValue( filename );*/

	std::ifstream file(filename);
	if (!file.good()) return false;

	std::string str((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	Document d;
	d.Parse(str.c_str());

	for (Value::ConstMemberIterator itr = d.MemberBegin(); itr != d.MemberEnd(); ++itr)         {             
		//printf("Member [%s] - type is [%s]\n", itr->name.GetString(), kTypeNames[itr->value.GetType()]);
		if(itr->value.IsObject()==true)
		{
			std::vector<unsigned int> aux;
			//std::cout<< "Process Object [" << itr->name.GetString() << "]" << std::endl;
			std::string s = itr->name.GetString();
			Value &v = d[s.c_str()];  // takes a "C" char*
			//if(v.IsObject())  std::cout<< "Its an object"<< std::endl;
			//if(v.IsArray())  std::cout<< "Its an Array"<< std::endl;
			//if(v.IsString())  std::cout<< "Its a String"<< std::endl;
			//if(v.IsNumber())  std::cout<< "Its a Number"<< std::endl;
			for (Value::ConstMemberIterator it = v.MemberBegin(); it != v.MemberEnd(); ++it)                 
			{  
				aux.push_back(it->value.GetInt());
			}
			m_muscleMatching[s] = aux;
		}
		//std::cout<<"\n"<< std::endl;
	}


	return true;
}

template <class DataTypes, class MassType>
bool MuscleActivation<DataTypes, MassType>::loadForceNodeFile(const char* filename)
{
	/*std::string meshFilename(filename);
	if (!sofa::helper::system::DataRepository.findFile (meshFilename))
	{
		serr << "Mesh \""<< filename <<"\" not found"<< sendl;
		return false;
	}
	this->f_fileNodes.setValue( filename );*/

	std::ifstream file(filename);
	if (!file.good()) return false;

	std::string line;
				

	while(std::getline(file, line)){
		int aux = std::stoi(line);
		m_forceNodes.push_back(aux);
	}

	return true;
}

template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::reinit()
{
    if (_topology && ((m_massDensity.getValue() > 0 && (vertexMassInfo.getValue().size() == 0 || edgeMassInfo.getValue().size() == 0)) || (m_massDensity.getValue()!= savedMass) ))
    {
        // resize array
        clear();

        /// prepare to store info in the vertex array
        vector<MassType>& my_vertexMassInfo = *vertexMassInfo.beginEdit();
        vector<MassType>& my_edgeMassInfo = *edgeMassInfo.beginEdit();

        unsigned int ndof = this->mstate->getSize();
        unsigned int nbEdges=_topology->getNbEdges();
        const helper::vector<topology::Edge>& edges = _topology->getEdges();

        my_vertexMassInfo.resize(ndof);
        my_edgeMassInfo.resize(nbEdges);

        const helper::vector< unsigned int > emptyAncestor;
        const helper::vector< double > emptyCoefficient;
        const helper::vector< helper::vector< unsigned int > > emptyAncestors;
        const helper::vector< helper::vector< double > > emptyCoefficients;

        // set vertex tensor to 0
        for (unsigned int i = 0; i<ndof; ++i)
            vertexMassHandler->applyCreateFunction(i, my_vertexMassInfo[i], emptyAncestor, emptyCoefficient);

        // set edge tensor to 0
        for (unsigned int i = 0; i<nbEdges; ++i)
            edgeMassHandler->applyCreateFunction(i, my_edgeMassInfo[i], edges[i], emptyAncestor, emptyCoefficient);

        // Create mass matrix depending on current Topology:
        if (_topology->getNbHexahedra()>0 && hexaGeo)  // Hexahedron topology
        {
            // create vector tensor by calling the hexahedron creation function on the entire mesh
            sofa::helper::vector<unsigned int> hexahedraAdded;
            setMassTopologyType(TOPOLOGY_HEXAHEDRONSET);
            int n = _topology->getNbHexahedra();
            for (int i = 0; i<n; ++i)
                hexahedraAdded.push_back(i);

            vertexMassHandler->applyHexahedronCreation(hexahedraAdded, _topology->getHexahedra(), emptyAncestors, emptyCoefficients);
            edgeMassHandler->applyHexahedronCreation(hexahedraAdded, _topology->getHexahedra(), emptyAncestors, emptyCoefficients);
            massLumpingCoeff = 2.5;
        }
        

		else if (_topology->getNbTetrahedra()>0 && bezierTetraGeo)  // Bezier Tetrahedron topology
        {
			vector<MassVector>& my_tetrahedronMassInfo = *tetrahedronMassInfo.beginEdit();


			size_t  nbTetrahedra=_topology->getNbTetrahedra();
			const helper::vector<topology::Tetra>& tetrahedra = _topology->getTetrahedra();

			my_tetrahedronMassInfo.resize(nbTetrahedra);
			 setMassTopologyType(TOPOLOGY_BEZIERTETRAHEDRONSET);
			// set vertex tensor to 0
			for (unsigned int i = 0; i<nbTetrahedra; ++i)
				tetrahedronMassHandler->applyCreateFunction(i, my_tetrahedronMassInfo[i], tetrahedra[i],emptyAncestor, emptyCoefficient);

            // create vector tensor by calling the tetrahedron creation function on the entire mesh
            sofa::helper::vector<unsigned int> tetrahedraAdded;
           

			size_t n = _topology->getNbTetrahedra();
			for (size_t i = 0; i<n; ++i)
				tetrahedraAdded.push_back(i);

//			tetrahedronMassHandler->applyTetrahedronCreation(tetrahedraAdded, _topology->getTetrahedra(), _topology->getTetrahedron(i),emptyAncestors, emptyCoefficients);
			massLumpingCoeff = 1.0;

			tetrahedronMassInfo.registerTopologicalData();
			tetrahedronMassInfo.endEdit();
        }
		else if (_topology->getNbTetrahedra()>0 && tetraGeo)  // Tetrahedron topology
        {
            // create vector tensor by calling the tetrahedron creation function on the entire mesh
            sofa::helper::vector<unsigned int> tetrahedraAdded;
            setMassTopologyType(TOPOLOGY_TETRAHEDRONSET);

            int n = _topology->getNbTetrahedra();
            for (int i = 0; i<n; ++i)
                tetrahedraAdded.push_back(i);

            vertexMassHandler->applyTetrahedronCreation(tetrahedraAdded, _topology->getTetrahedra(), emptyAncestors, emptyCoefficients);
            edgeMassHandler->applyTetrahedronCreation(tetrahedraAdded, _topology->getTetrahedra(), emptyAncestors, emptyCoefficients);
            massLumpingCoeff = 2.5;
        }
        else if (_topology->getNbQuads()>0 && quadGeo)  // Quad topology
        {
            // create vector tensor by calling the quad creation function on the entire mesh
            sofa::helper::vector<unsigned int> quadsAdded;
            setMassTopologyType(TOPOLOGY_QUADSET);

            int n = _topology->getNbQuads();
            for (int i = 0; i<n; ++i)
                quadsAdded.push_back(i);

            vertexMassHandler->applyQuadCreation(quadsAdded, _topology->getQuads(), emptyAncestors, emptyCoefficients);
            edgeMassHandler->applyQuadCreation(quadsAdded, _topology->getQuads(), emptyAncestors, emptyCoefficients);
            massLumpingCoeff = 2.0;
        }
        else if (_topology->getNbTriangles()>0 && triangleGeo) // Triangle topology
        {
            // create vector tensor by calling the triangle creation function on the entire mesh
            sofa::helper::vector<unsigned int> trianglesAdded;
            setMassTopologyType(TOPOLOGY_TRIANGLESET);

            int n = _topology->getNbTriangles();
            for (int i = 0; i<n; ++i)
                trianglesAdded.push_back(i);

            vertexMassHandler->applyTriangleCreation(trianglesAdded, _topology->getTriangles(), emptyAncestors, emptyCoefficients);
            edgeMassHandler->applyTriangleCreation(trianglesAdded, _topology->getTriangles(), emptyAncestors, emptyCoefficients);
            massLumpingCoeff = 2.0;
        }

        vertexMassInfo.registerTopologicalData();
        edgeMassInfo.registerTopologicalData();

        vertexMassInfo.endEdit();
        edgeMassInfo.endEdit();
    }
}


template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::copyVertexMass() {}


template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::clear()
{
    MassVector& vertexMass = *vertexMassInfo.beginEdit();
    MassVector& edgeMass = *edgeMassInfo.beginEdit();
	MassVectorVector& tetrahedronMass = *tetrahedronMassInfo.beginEdit();
    vertexMass.clear();
    edgeMass.clear();
	tetrahedronMass.clear();
    vertexMassInfo.endEdit();
    edgeMassInfo.endEdit();
	tetrahedronMassInfo.endEdit();
}


// -- Mass interface
template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::addMDx(const core::MechanicalParams*, DataVecDeriv& vres, const DataVecDeriv& vdx, SReal factor)
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    const MassVector &edgeMass= edgeMassInfo.getValue();

    helper::WriteAccessor< DataVecDeriv > res = vres;
    helper::ReadAccessor< DataVecDeriv > dx = vdx;

    SReal massTotal = 0.0;

    //using a lumped matrix (default)-----
    if(this->lumping.getValue())
    {
        for (size_t i=0; i<dx.size(); i++)
        {
            res[i] += dx[i] * vertexMass[i] * massLumpingCoeff * (Real)factor;
            massTotal += vertexMass[i]*massLumpingCoeff * (Real)factor;
        }
        
    }


    //using a sparse matrix---------------
	else if (getMassTopologyType()!=TOPOLOGY_BEZIERTETRAHEDRONSET) 
	{
		size_t nbEdges=_topology->getNbEdges();
		size_t v0,v1;

		for (unsigned int i=0; i<dx.size(); i++)
		{
			res[i] += dx[i] * vertexMass[i] * (Real)factor;
			massTotal += vertexMass[i] * (Real)factor;
		}

		Real tempMass=0.0;

		for (unsigned int j=0; j<nbEdges; ++j)
		{
			tempMass = edgeMass[j] * (Real)factor;

			v0=_topology->getEdge(j)[0];
			v1=_topology->getEdge(j)[1];

			res[v0] += dx[v1] * tempMass;
			res[v1] += dx[v0] * tempMass;

			massTotal += 2*edgeMass[j] * (Real)factor;
		}
	} else if (bezierTetraGeo ){
			topology::BezierDegreeType degree=bezierTetraGeo->getTopologyContainer()->getDegree();
			size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
			topology::VecPointID indexArray;
			size_t nbTetras=_topology->getNbTetrahedra();
#ifdef NDEBUG
			assert(tetrahedronMassInfo.size()==(nbControlPoints*(nbControlPoints+1)/2));
#endif
			// go through the mass stored in each tetrahedron element
			size_t rank=0;
			MassType tempMass;
			size_t v0,v1;
			// loop over each tetrahedron of size nbControlPoints*nbControlPoints
			for (size_t i=0; i<nbTetras; i++) {
				indexArray.clear();
				/// get the global index of each control point in the tetrahedron
				bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfBezierPointsInTetrahedron(i,indexArray) ;
				// get the mass matrix in the tetrahedron
				const MassVector &mv=tetrahedronMassInfo.getValue()[i];
				// loop over each entry in the mass matrix of size nbControlPoints*(nbControlPoints+1)/2
				rank=0;
				for (size_t j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for (size_t k=j; k<nbControlPoints; ++k,++rank) {
						v1 = indexArray[k];
						tempMass =mv[rank] * (Real)factor;					
						if (k>j) {
							res[v0] += dx[v1] * tempMass;
							res[v1] += dx[v0] * tempMass;
							massTotal += 2*tempMass;
						} else {
							res[v0] += dx[v0] * tempMass;
							massTotal += tempMass;
						}
					}
				}
			}
		}
	if(printMass.getValue() && (this->getContext()->getTime()==0.0))
        sout<<"Total Mass = "<<massTotal<<sendl;

	if(printMass.getValue())
	{
		std::map < std::string, sofa::helper::vector<double> >& graph = *f_graph.beginEdit();
		sofa::helper::vector<double>& graph_error = graph["Mass variations"];
		graph_error.push_back(massTotal+0.000001);

		f_graph.endEdit();
	}

        
    

}



template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::accFromF(const core::MechanicalParams*, DataVecDeriv& a, const DataVecDeriv& f)
{
    helper::WriteAccessor< DataVecDeriv > _a = a;
    const VecDeriv& _f = f.getValue();
    const MassVector &vertexMass= vertexMassInfo.getValue();

    if(this->lumping.getValue())
    {
        for (unsigned int i=0; i<vertexMass.size(); i++)
        {
            _a[i] = _f[i] / ( vertexMass[i] * massLumpingCoeff);
        }
    }
    else
    {
        (void)a;
        (void)f;
        serr << "WARNING: the methode 'accFromF' can't be used with MuscleActivation as this SPARSE mass matrix can't be inversed easily. \nPlease proceed to mass lumping." << sendl;
        return;
    }
}




#ifdef SOFA_SUPPORT_MOVING_FRAMES
template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::addForce(const core::MechanicalParams*, DataVecDeriv& vf, const DataVecCoord& vx, const DataVecDeriv& vv)
{
    helper::WriteAccessor< DataVecDeriv > f = vf;
    helper::ReadAccessor< DataVecCoord > x = vx;
    helper::ReadAccessor< DataVecDeriv > v = vv;

    const MassVector &vertexMass= vertexMassInfo.getValue();

    // gravity
    Vec3d g ( this->getContext()->getGravity() );
    Deriv theGravity;
    DataTypes::set ( theGravity, g[0], g[1], g[2]);

    // velocity-based stuff
    core::objectmodel::BaseContext::SpatialVector vframe = this->getContext()->getVelocityInWorld();
    core::objectmodel::BaseContext::Vec3 aframe = this->getContext()->getVelocityBasedLinearAccelerationInWorld() ;

    // project back to local frame
    vframe = this->getContext()->getPositionInWorld() / vframe;
    aframe = this->getContext()->getPositionInWorld().backProjectVector( aframe );

    // add weight and inertia force
    if(this->m_separateGravity.getValue())
        for (unsigned int i=0; i<x.size(); ++i)
            f[i] += massLumpingCoeff + core::behavior::inertiaForce(vframe,aframe,vertexMass[i] * massLumpingCoeff ,x[i],v[i]);
    else for (unsigned int i=0; i<x.size(); ++i)
            f[i] += theGravity * vertexMass[i] * massLumpingCoeff + core::behavior::inertiaForce(vframe,aframe,vertexMass[i] * massLumpingCoeff ,x[i],v[i]);
}
#else
template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::addForce(const core::MechanicalParams* /*mp*/, DataVecDeriv& vf, const DataVecCoord& d_x, const DataVecDeriv& )
{
	helper::WriteAccessor< DataVecDeriv > f = vf;
	const VecCoord& x = d_x.getValue();

	//std::cout << "f1: " << f[0] << " " << f[1] << " " << f[2] << std::endl;
    //if gravity was added separately (in solver's "solve" method), then nothing to do here
    if(this->m_separateGravity.getValue())
        return ;

    const MassVector &vertexMass= vertexMassInfo.getValue();

    // gravity
    defaulttype::Vec3d g ( this->getContext()->getGravity() );
    Deriv theGravity;
    DataTypes::set ( theGravity, g[0], g[1], g[2]);
	
    // add weight and inertia force
    for (unsigned int i=0; i<f.size(); ++i){
		if(i  < 3){
			//std::cout << "a: " << theGravity * vertexMass[i] * massLumpingCoeff << " " << theGravity << " " << vertexMass[i] << " " << massLumpingCoeff << " ";
		}
        f[i] += theGravity * vertexMass[i] * massLumpingCoeff;
	}


	counter = this->getContext()->getTime();
	
	int numVert = _topology->getNbPoints();

	for(unsigned int i = 0; i < m_forceHelper.size(); i++){
		m_forceHelper[i](0) = 0;
		m_forceHelper[i](1) = 0;
		m_forceHelper[i](2) = 0;
	}
	

	//std::cout << "mp: " << mechanicalObject->rotation.getValue()[0] << std::endl;
	
	//                      | Y
	//              1       |        2
	//                      |
	//          0           |            3       Z
	//---------------------------------------------
	//          7           |            4 
	//                      | 
	//              6       |        5*
	//                      |

	//* row with 11 muscles

	if (!f_signalFile.getValue().empty())
	{
		Real neuronTimestep = f_neuronTimestep.getValue();
		unsigned int timestepNum = floor(this->getContext()->getDt()/neuronTimestep);
		unsigned int startingTimestep = floor(this->getContext()->getTime()/neuronTimestep);
		const double pi=3.14159265358979323846264338327950288;

		if(m_activationSignals.size() <= startingTimestep){
			return;
		}

		std::vector<std::vector<double>> currentSignals;
		for(unsigned int i = 0; i < m_orderedMuscles.size(); i++){
			std::vector<double> aux;
			for(unsigned int j = 0; j < m_orderedMuscles.at(i).size(); j++){
				aux.push_back(0.0);
			}
			currentSignals.push_back(aux);
		}
	
		for(unsigned int i = 0; i < timestepNum; i++){
			unsigned int j = 0;
			for(std::map<std::string,std::vector<unsigned int>>::iterator it = m_muscleMatching.begin(); it != m_muscleMatching.end(); it++){
				//take the last signal. Insert sum, interpolation or whatever here. (with the last timesteps)

				//taking the last
				/*if(m_activationSignals[startingTimestep + i][j] == 0){
					currentSignals.at(it->second.at(0)).at(it->second.at(1)) = 0;
				}
				else{
					currentSignals.at(it->second.at(0)).at(it->second.at(1)) = m_activationSignals[startingTimestep + i][j];
				}*/

				//average value
				currentSignals.at(it->second.at(0)).at(it->second.at(1)) += m_activationSignals[startingTimestep + i][j];
				j++;
			}
		}

		//average value
		for(unsigned int i = 0; i < currentSignals.size(); i++){
			for(unsigned int j = 0; j < currentSignals.at(i).size(); j++){
				currentSignals.at(i).at(j) = currentSignals.at(i).at(j)/timestepNum;
			}
		}

		for(unsigned int i = 0; i < m_orderedMuscles.size(); i++){
			for(unsigned int j = 0; j < m_orderedMuscles.at(i).size(); j++){
				//if(j > 0) continue;
				double muscleMidRef = (m_orderedMuscles.at(i).at(j).size() - 1)/2;
				for(unsigned int k = 0; k < m_orderedMuscles.at(i).at(j).size(); k++){
					if(k < muscleMidRef){
						for(unsigned int l = 0; l < m_orderedMuscles.at(i).at(j).at(k).size(); l++){
							if(l == 0 && k == 0){
								m_saveSignals.at(i).at(j) = currentSignals.at(i).at(j)*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							}
							f[m_orderedMuscles[i][j][k][l]][0] += currentSignals.at(i).at(j)*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += currentSignals.at(i).at(j)*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += currentSignals.at(i).at(j)*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							//if(m_orderedMuscles[i][j][k][l] == 331 || m_orderedMuscles[i][j][k][l] == 348) std::cout << "f2: " << counter << " ; " << f[331] << " , " << f[348] << std::endl;
						}																			  
					}																				  
					else if(k > muscleMidRef){														  
						for(unsigned int l = 0; l < m_orderedMuscles.at(i).at(j).at(k).size(); l++){  
							f[m_orderedMuscles[i][j][k][l]][0] += currentSignals.at(i).at(j)*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += currentSignals.at(i).at(j)*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += currentSignals.at(i).at(j)*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							//if(m_orderedMuscles[i][j][k][l] == 331 || m_orderedMuscles[i][j][k][l] == 348) std::cout << "f3: " << counter << " ; " << f[331] << " , " << f[348] << std::endl;
						}
					}	
				}
			}
		}
	}
	else{
		// Centro de masas
		Eigen::Vector3d cm;
		if(f_correctiveForces.getValue()){
			cm = centerMass(x, numVert);
			//std::cout << "cm: " << cm << std::endl;
		}

		Eigen::Vector3d sumForce; // sumForce = sum(Finterna_i)
		sumForce(0) = sumForce(1) = sumForce(2) = 0.0;
		Eigen::Vector3d sumTorque; // sumTorque = sum(torque_i)
		sumTorque(0) = sumTorque(1) = sumTorque(2) = 0.0;
		Eigen::Vector3d pos;
		Eigen::Vector3d lambdaForce, lambdaTorque, position, forc, tor;
		const double pi=3.14159265358979323846264338327950288;
		std::vector<int> sign;
		sign.push_back(1);
		sign.push_back(1);
		sign.push_back(-1);
		sign.push_back(-1);
		sign.push_back(-1);
		sign.push_back(-1);
		sign.push_back(1);
		sign.push_back(1);

		
		double direction = 1.0;
		//check if we are in the first half of the period or not
		double check = counter - floor(counter/(m_activationPeriod1+m_activationPeriod2+m_activationStop1+m_activationStop2))*(m_activationPeriod1+m_activationPeriod2+m_activationStop1+m_activationStop2);
		for(unsigned int i = 0; i < m_orderedMuscles.size(); i++){
			/*if(counter < 460){
				break;
			}*/
			for(unsigned int j = m_minActivatedMuscles; j < m_orderedMuscles.at(i).size() && j < m_maxActivatedMuscles; j++){
				double muscleMidRef = (m_orderedMuscles.at(i).at(j).size() - 1)/2;
				if(check < m_activationPeriod1){
					m_saveSignals.at(i).at(j) = sign[i]*sin(check*pi*2/m_activationPeriod1)*sin(3*pi*j/24);
				}
				else if(check < m_activationPeriod1+m_activationStop1){
					m_saveSignals.at(i).at(j) = 0;
				}
				else if(check < m_activationPeriod1+m_activationStop1+m_activationPeriod2){
					m_saveSignals.at(i).at(j) = -sign[i]*sin((check-m_activationPeriod1-m_activationStop1)*pi*2/m_activationPeriod2)*sin(3*pi*j/24);
				}
				else if(check < m_activationPeriod1+m_activationStop1+m_activationPeriod2+m_activationStop2){
					m_saveSignals.at(i).at(j) = 0;
				}
				for(unsigned int k = 0; k < m_orderedMuscles.at(i).at(j).size(); k++){
					if(k < muscleMidRef){
						for(unsigned int l = 0; l < m_orderedMuscles.at(i).at(j).at(k).size(); l++){
							//everything
							/*f[m_orderedMuscles[i][j][k][l]][0] += sign[i]*sin(counter*pi*2/m_activationPeriod - 3*pi*j/24)*(muscleMidRef-k)*activationConstant*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += sign[i]*sin(counter*pi*2/m_activationPeriod - 3*pi*j/24)*(muscleMidRef-k)*activationConstant*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += sign[i]*sin(counter*pi*2/m_activationPeriod - 3*pi*j/24)*(muscleMidRef-k)*activationConstant*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);*/

							//all the row the same
							//f[m_orderedMuscles[i][j][k][l]][0] += sign[i]*sin(counter*pi*2/m_activationPeriod)*(muscleMidRef-k)*m_activationConstant*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							//f[m_orderedMuscles[i][j][k][l]][1] += sign[i]*sin(counter*pi*2/m_activationPeriod)*(muscleMidRef-k)*m_activationConstant*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							//f[m_orderedMuscles[i][j][k][l]][2] += sign[i]*sin(counter*pi*2/m_activationPeriod)*(muscleMidRef-k)*m_activationConstant*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							//if(m_orderedMuscles[i][j][k][l] == 331 || m_orderedMuscles[i][j][k][l] == 348) std::cout << "ei: " << x.size() << " " << m_orderedMuscles[i][j][k+1][l] << " " << x[m_orderedMuscles[i][j][k+1][l]][0] << " " << m_orderedMuscles[i][j][k][l] << " " << x[m_orderedMuscles[i][j][k][l]][0] << " " << x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0] << " " << x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1] << " " << x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2] << std::endl;
							//if(m_orderedMuscles[i][j][k][l] == 331 || m_orderedMuscles[i][j][k][l] == 348) std::cout << "f2: " << counter << " ; " <<  f[331] << " , " << f[348] << std::endl;
							//if(i == 0 && j==m_minActivatedMuscles && k==0 && l==0) std::cout << "check: " << check << "\n";
							
							/*if(check < m_activationPeriod1){
								//if(i == 0 && j==m_minActivatedMuscles && k==0 && l==0) std::cout << "period 1: " << std::endl;
								f[m_orderedMuscles[i][j][k][l]][0] += sign[i]*sin(check*pi/m_activationPeriod1)*sin(3*pi*j/24)*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
								f[m_orderedMuscles[i][j][k][l]][1] += sign[i]*sin(check*pi/m_activationPeriod1)*sin(3*pi*j/24)*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
								f[m_orderedMuscles[i][j][k][l]][2] += sign[i]*sin(check*pi/m_activationPeriod1)*sin(3*pi*j/24)*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							}
							else if(check < m_activationPeriod1+m_activationStop1){
								//if(i == 0 && j==m_minActivatedMuscles && k==0 && l==0) std::cout << "stop 1: " << std::endl;
								f[m_orderedMuscles[i][j][k][l]][0] += 0;
								f[m_orderedMuscles[i][j][k][l]][1] += 0;
								f[m_orderedMuscles[i][j][k][l]][2] += 0;
							}
							else if(check < m_activationPeriod1+m_activationStop1+m_activationPeriod2){
								//if(i == 0 && j==m_minActivatedMuscles && k==0 && l==0) std::cout << "period 2: " << std::endl;
								f[m_orderedMuscles[i][j][k][l]][0] -= sign[i]*sin((check-m_activationPeriod1-m_activationStop1)*pi/m_activationPeriod2)*sin(3*pi*j/24)*(muscleMidRef-k)*m_activationConstant2*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
								f[m_orderedMuscles[i][j][k][l]][1] -= sign[i]*sin((check-m_activationPeriod1-m_activationStop1)*pi/m_activationPeriod2)*sin(3*pi*j/24)*(muscleMidRef-k)*m_activationConstant2*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
								f[m_orderedMuscles[i][j][k][l]][2] -= sign[i]*sin((check-m_activationPeriod1-m_activationStop1)*pi/m_activationPeriod2)*sin(3*pi*j/24)*(muscleMidRef-k)*m_activationConstant2*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							}
							else if(check < m_activationPeriod1+m_activationStop1+m_activationPeriod2+m_activationStop2){
								//if(i == 0 && j==m_minActivatedMuscles && k==0 && l==0) std::cout << "stop 2: " << std::endl;
								f[m_orderedMuscles[i][j][k][l]][0] += 0;
								f[m_orderedMuscles[i][j][k][l]][1] += 0;
								f[m_orderedMuscles[i][j][k][l]][2] += 0;
							}*/

							//if(i == 0 && j==0 && k == 0 && l==0)std::cout << "ei: " << (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction)+1.f)/2 ) << std::endl;
							/*f[m_orderedMuscles[i][j][k][l]][0] += (float)( (sin(3*pi*j/24 - sign[i]*pi/2 - 2*pi*counter*direction/m_activationPeriod1)+1.f)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += (float)( (sin(3*pi*j/24 - sign[i]*pi/2 - 2*pi*counter*direction/m_activationPeriod1)+1.f)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += (float)( (sin(3*pi*j/24 - sign[i]*pi/2 - 2*pi*counter*direction/m_activationPeriod1)+1.f)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);*/
							
							/*f[m_orderedMuscles[i][j][k][l]][0] += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction) - 1)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction) - 1)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction) - 1)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);*/
							
							m_forceHelper[m_orderedMuscles[i][j][k][l]](0) += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction)+1.f)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							m_forceHelper[m_orderedMuscles[i][j][k][l]](1) += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction)+1.f)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							m_forceHelper[m_orderedMuscles[i][j][k][l]](2) += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction)+1.f)/2 )*(muscleMidRef-k)*m_activationConstant1*(x[m_orderedMuscles[i][j][k+1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							if(f_correctiveForces.getValue()){
								pos(0) = x[m_orderedMuscles[i][j][k][l]][0]; pos(1) = x[m_orderedMuscles[i][j][k][l]][1]; pos(2) = x[m_orderedMuscles[i][j][k][l]][2];
								sumForce += m_forceHelper[m_orderedMuscles[i][j][k][l]]; //F0
								sumTorque += (pos - cm).cross(m_forceHelper[m_orderedMuscles[i][j][k][l]]); //T0
							}
						}																			  
					}																				  
					else if(k > muscleMidRef){														  
						for(unsigned int l = 0; l < m_orderedMuscles.at(i).at(j).at(k).size(); l++){
							//everything
							/*f[m_orderedMuscles[i][j][k][l]][0] += sign[i]*sin(counter*pi*2/m_activationPeriod - 3*pi*j/24)*(k-muscleMidRef)*activationConstant*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += sign[i]*sin(counter*pi*2/m_activationPeriod - 3*pi*j/24)*(k-muscleMidRef)*activationConstant*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += sign[i]*sin(counter*pi*2/m_activationPeriod - 3*pi*j/24)*(k-muscleMidRef)*activationConstant*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);*/

							//all the row the same
							//f[m_orderedMuscles[i][j][k][l]][0] += sign[i]*sin(counter*pi*2/m_activationPeriod)*(k-muscleMidRef)*m_activationConstant*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							//f[m_orderedMuscles[i][j][k][l]][1] += sign[i]*sin(counter*pi*2/m_activationPeriod)*(k-muscleMidRef)*m_activationConstant*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							//f[m_orderedMuscles[i][j][k][l]][2] += sign[i]*sin(counter*pi*2/m_activationPeriod)*(k-muscleMidRef)*m_activationConstant*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							//if(m_orderedMuscles[i][j][k][l] == 331 || m_orderedMuscles[i][j][k][l] == 348) std::cout << "f3: " << counter << " ; " << f[331] << " , " << f[348] << std::endl;

							/*if(check < m_activationPeriod1){
								f[m_orderedMuscles[i][j][k][l]][0] += sign[i]*sin(check*pi/m_activationPeriod1)*sin(3*pi*j/24)*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
								f[m_orderedMuscles[i][j][k][l]][1] += sign[i]*sin(check*pi/m_activationPeriod1)*sin(3*pi*j/24)*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
								f[m_orderedMuscles[i][j][k][l]][2] += sign[i]*sin(check*pi/m_activationPeriod1)*sin(3*pi*j/24)*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							}
							else if(check < m_activationPeriod1+m_activationStop1){
								f[m_orderedMuscles[i][j][k][l]][0] += 0;
								f[m_orderedMuscles[i][j][k][l]][1] += 0;
								f[m_orderedMuscles[i][j][k][l]][2] += 0;
							}
							else if(check < m_activationPeriod1+m_activationStop1+m_activationPeriod2){
								f[m_orderedMuscles[i][j][k][l]][0] -= sign[i]*sin((check-m_activationPeriod1)*pi/m_activationPeriod2)*sin(3*pi*j/24)*(k-muscleMidRef)*m_activationConstant2*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
								f[m_orderedMuscles[i][j][k][l]][1] -= sign[i]*sin((check-m_activationPeriod1)*pi/m_activationPeriod2)*sin(3*pi*j/24)*(k-muscleMidRef)*m_activationConstant2*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
								f[m_orderedMuscles[i][j][k][l]][2] -= sign[i]*sin((check-m_activationPeriod1)*pi/m_activationPeriod2)*sin(3*pi*j/24)*(k-muscleMidRef)*m_activationConstant2*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);
							}
							else if(check < m_activationPeriod1+m_activationStop1+m_activationPeriod2+m_activationStop2){
								f[m_orderedMuscles[i][j][k][l]][0] += 0;
								f[m_orderedMuscles[i][j][k][l]][1] += 0;
								f[m_orderedMuscles[i][j][k][l]][2] += 0;
							}*/

							/*f[m_orderedMuscles[i][j][k][l]][0] += (float)( (sin(3*pi*j/24 - sign[i]*pi/2 - 2*pi*counter*direction/m_activationPeriod1)+1.f)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += (float)( (sin(3*pi*j/24 - sign[i]*pi/2 - 2*pi*counter*direction/m_activationPeriod1)+1.f)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += (float)( (sin(3*pi*j/24 - sign[i]*pi/2 - 2*pi*counter*direction/m_activationPeriod1)+1.f)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);*/

							/*f[m_orderedMuscles[i][j][k][l]][0] += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction) - 1)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							f[m_orderedMuscles[i][j][k][l]][1] += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction) - 1)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							f[m_orderedMuscles[i][j][k][l]][2] += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction) - 1)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);*/

							m_forceHelper[m_orderedMuscles[i][j][k][l]](0) += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction)+1.f)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][0]-x[m_orderedMuscles[i][j][k][l]][0]);
							m_forceHelper[m_orderedMuscles[i][j][k][l]](1) += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction)+1.f)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][1]-x[m_orderedMuscles[i][j][k][l]][1]);
							m_forceHelper[m_orderedMuscles[i][j][k][l]](2) += (float)( 0.4f * (sin(3*pi*j/24 + sign[i]*pi/2 + pi/2 - 0.5*counter*direction)+1.f)/2 )*(k-muscleMidRef)*m_activationConstant1*(x[m_orderedMuscles[i][j][k-1][l]][2]-x[m_orderedMuscles[i][j][k][l]][2]);

							if(f_correctiveForces.getValue()){
								pos(0) = x[m_orderedMuscles[i][j][k][l]][0]; pos(1) = x[m_orderedMuscles[i][j][k][l]][1]; pos(2) = x[m_orderedMuscles[i][j][k][l]][2];
								sumForce += m_forceHelper[m_orderedMuscles[i][j][k][l]]; //F0
								sumTorque += (pos - cm).cross(m_forceHelper[m_orderedMuscles[i][j][k][l]]); //T0
							}
						}
					}
				}
			}
		}
		// Calculo de lambdas
		if(f_correctiveForces.getValue()){
			multLagrange(lambdaForce, lambdaTorque, x, sumForce, sumTorque, cm, numVert);
		}
		tor(0) = tor(1) = tor(2) = 0.0;
		forc(0) = forc(1) = forc(2) = 0.0;

		for(unsigned int i = 0; i < m_orderedMuscles.size(); i++){
			for(unsigned int j = m_minActivatedMuscles; j < m_orderedMuscles.at(i).size() && j < m_maxActivatedMuscles; j++){
				for(unsigned int k = 0; k < m_orderedMuscles.at(i).at(j).size(); k++){
					for(unsigned int l = 0; l < m_orderedMuscles.at(i).at(j).at(k).size(); l++){
						if(f_correctiveForces.getValue()){

							// Calculo de fuerzas correctivas

							position(0) = x[m_orderedMuscles[i][j][k][l]][0];		position(1) = x[m_orderedMuscles[i][j][k][l]][1];		position(2) = x[m_orderedMuscles[i][j][k][l]][2];

							Eigen::Vector3d F =	createForceCorrective(position, cm, lambdaForce,  lambdaTorque);

							f[m_orderedMuscles[i][j][k][l]][0] += m_forceHelper[m_orderedMuscles[i][j][k][l]][0] + F(0);
							f[m_orderedMuscles[i][j][k][l]][1] += m_forceHelper[m_orderedMuscles[i][j][k][l]][1] + F(1);
							f[m_orderedMuscles[i][j][k][l]][2] += m_forceHelper[m_orderedMuscles[i][j][k][l]][2] + F(2);

							forc += F;
							tor += (position - cm).cross(F);
						}
						else{
							f[m_orderedMuscles[i][j][k][l]][0] += m_forceHelper[m_orderedMuscles[i][j][k][l]][0];
							f[m_orderedMuscles[i][j][k][l]][1] += m_forceHelper[m_orderedMuscles[i][j][k][l]][1];
							f[m_orderedMuscles[i][j][k][l]][2] += m_forceHelper[m_orderedMuscles[i][j][k][l]][2];
						}
					}
				}
			}
		}
	}
	m_saveSignalsOrdered.clear();
	for(std::map<std::string,std::vector<unsigned int>>::iterator it = m_muscleMatching.begin(); it != m_muscleMatching.end(); it++){
		m_saveSignalsOrdered.push_back(m_saveSignals.at(it->second.at(0)).at(it->second.at(1)));
	}
	//std::cout << "f3: " << f[0] << " " << f[1] << " " << f[2] << std::endl;
}
#endif


template <class DataTypes, class MassType>
SReal MuscleActivation<DataTypes, MassType>::getKineticEnergy( const core::MechanicalParams*, const DataVecDeriv& vv ) const
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    const MassVector &edgeMass= edgeMassInfo.getValue();

    helper::ReadAccessor< DataVecDeriv > v = vv;

    unsigned int nbEdges=_topology->getNbEdges();
    unsigned int v0,v1;

    SReal e = 0;

    for (unsigned int i=0; i<v.size(); i++)
    {
        e += dot(v[i],v[i]) * vertexMass[i]; // v[i]*v[i]*masses[i] would be more efficient but less generic
    }
	if (getMassTopologyType()!=TOPOLOGY_BEZIERTETRAHEDRONSET) {
		for (unsigned int i=0; i<nbEdges; ++i)
		{
			v0=_topology->getEdge(i)[0];
			v1=_topology->getEdge(i)[1];

			e += 2*dot(v[v0],v[v1])*edgeMass[i];

		} 
	} else if (bezierTetraGeo ){
			topology::BezierDegreeType degree=bezierTetraGeo->getTopologyContainer()->getDegree();
			size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
			topology::VecPointID indexArray;
			size_t nbTetras=_topology->getNbTetrahedra();
#ifdef NDEBUG
			assert(tetrahedronMassInfo.size()==(nbControlPoints*(nbControlPoints+1)/2));
#endif
			// go through the mass stored in each tetrahedron element
			size_t rank=0;
			// loop over each tetrahedron of size nbControlPoints*nbControlPoints
			for (size_t i=0; i<nbTetras; i++) {
				indexArray.clear();
				/// get the global index of each control point in the tetrahedron
				bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfBezierPointsInTetrahedron(i,indexArray) ;
				// get the mass matrix in the tetrahedron
				const MassVector &mv=tetrahedronMassInfo.getValue()[i];
			//	MassVector mv;
				// loop over each entry in the mass matrix of size nbControlPoints*(nbControlPoints+1)/2
				for (size_t j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for (size_t k=j; k<nbControlPoints; ++k,++rank) {
						v1 = indexArray[k];
						
						if (k>j) {
							e += 2*dot(v[v0],v[v1])*mv[rank];
						} else 
							e += dot(v[v0],v[v1])*mv[rank];
					}
				}
			}
		}


    return e/2;
}


template <class DataTypes, class MassType>
SReal MuscleActivation<DataTypes, MassType>::getPotentialEnergy( const core::MechanicalParams*, const DataVecCoord& vx) const
{
    const MassVector &vertexMass= vertexMassInfo.getValue();

    helper::ReadAccessor< DataVecCoord > x = vx;

    SReal e = 0;
    // gravity
    defaulttype::Vec3d g ( this->getContext()->getGravity() );
    Deriv theGravity;
    DataTypes::set ( theGravity, g[0], g[1], g[2]);

    for (unsigned int i=0; i<x.size(); i++)
        e -= dot(theGravity,x[i])*vertexMass[i] * massLumpingCoeff;

    return e;
}


// does nothing by default, need to be specialized in .cpp
template <class DataTypes, class MassType>
defaulttype::Vector6 MuscleActivation<DataTypes, MassType>::getMomentum ( const core::MechanicalParams*, const DataVecCoord& /*vx*/, const DataVecDeriv& /*vv*/  ) const
{
    return defaulttype::Vector6();
}



template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::addGravityToV(const core::MechanicalParams* mparams, DataVecDeriv& d_v)
{
    if(this->mstate && mparams)
    {
        VecDeriv& v = *d_v.beginEdit();

        // gravity
        defaulttype::Vec3d g ( this->getContext()->getGravity() );
        Deriv theGravity;
        DataTypes::set ( theGravity, g[0], g[1], g[2]);
        Deriv hg = theGravity * (typename DataTypes::Real)(mparams->dt());

        for (unsigned int i=0; i<v.size(); i++)
            v[i] += hg;
        d_v.endEdit();
    }

}


template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::addMToMatrix(const core::MechanicalParams *mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    const MassVector &edgeMass= edgeMassInfo.getValue();

    size_t nbEdges=_topology->getNbEdges();
    size_t v0,v1;

    const int N = defaulttype::DataTypeInfo<Deriv>::size();
    AddMToMatrixFunctor<Deriv,MassType> calc;
    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    sofa::defaulttype::BaseMatrix* mat = r.matrix;
    Real mFactor = (Real)mparams->mFactorIncludingRayleighDamping(this->rayleighMass.getValue());

    if((int)mat->colSize() != (_topology->getNbPoints()*N) || (int)mat->rowSize() != (_topology->getNbPoints()*N))
    {
        serr<<"Wrong size of the input Matrix: need resize in addMToMatrix function."<<sendl;
        mat->resize(_topology->getNbPoints()*N,_topology->getNbPoints()*N);
    }

    SReal massTotal=0.0;

    if(this->lumping.getValue())
    {
        for (size_t i=0; i<vertexMass.size(); i++)
        {
            calc(r.matrix, vertexMass[i] * massLumpingCoeff, r.offset + N*i, mFactor);
            massTotal += vertexMass[i] * massLumpingCoeff;
        }

        if(printMass.getValue() && (this->getContext()->getTime()==0.0))
            //std::cout<<"Total Mass = "<<massTotal<<std::endl;

        if(printMass.getValue())
        {
            std::map < std::string, sofa::helper::vector<double> >& graph = *f_graph.beginEdit();
            sofa::helper::vector<double>& graph_error = graph["Mass variations"];
            graph_error.push_back(massTotal);

            f_graph.endEdit();
        }
    }


    else
    {
		if (getMassTopologyType()!=TOPOLOGY_BEZIERTETRAHEDRONSET) {
			for (size_t i=0; i<vertexMass.size(); i++)
			{
				calc(r.matrix, vertexMass[i], r.offset + N*i, mFactor);
				massTotal += vertexMass[i];
			}


			for (size_t j=0; j<nbEdges; ++j)
			{
				v0=_topology->getEdge(j)[0];
				v1=_topology->getEdge(j)[1];

				calc(r.matrix, edgeMass[j], r.offset + N*v0, r.offset + N*v1, mFactor);
				calc(r.matrix, edgeMass[j], r.offset + N*v1, r.offset + N*v0, mFactor);

				massTotal += 2*edgeMass[j];
			}
		} else if (bezierTetraGeo ){
			topology::BezierDegreeType degree=bezierTetraGeo->getTopologyContainer()->getDegree();
			size_t nbControlPoints=(degree+1)*(degree+2)*(degree+3)/6;
			topology::VecPointID indexArray;
			size_t nbTetras=_topology->getNbTetrahedra();
#ifdef NDEBUG
			assert(tetrahedronMassInfo.size()==(nbControlPoints*(nbControlPoints+1)/2));
#endif
			// go through the mass stored in each tetrahedron element
			size_t rank=0;
			// loop over each tetrahedron of size nbControlPoints*nbControlPoints
			for (size_t i=0; i<nbTetras; i++) {
				indexArray.clear();
				/// get the global index of each control point in the tetrahedron
				bezierTetraGeo->getTopologyContainer()->getGlobalIndexArrayOfBezierPointsInTetrahedron(i,indexArray) ;
				// get the mass matrix in the tetrahedron
				MassVector &mv=tetrahedronMassInfo[i];
				// loop over each entry in the mass matrix of size nbControlPoints*(nbControlPoints+1)/2
				for (size_t j=0; j<nbControlPoints; ++j) {
					v0 = indexArray[j];
					for (size_t k=j; k<nbControlPoints; ++k,++rank) {
						v1 = indexArray[k];
						calc(r.matrix, mv[rank], r.offset + N*v0, r.offset + N*v1, mFactor);
						
						if (k>j) {
							calc(r.matrix, mv[rank], r.offset + N*v1, r.offset + N*v0, mFactor);
							massTotal += 2*mv[rank];
						} else 
							massTotal += mv[rank];
					}
				}
			}
		}
        if(printMass.getValue() && (this->getContext()->getTime()==0.0))
            //std::cout<<"Total Mass  = "<<massTotal<<std::endl;

        if(printMass.getValue())
        {
            std::map < std::string, sofa::helper::vector<double> >& graph = *f_graph.beginEdit();
            sofa::helper::vector<double>& graph_error = graph["Mass variations"];
            graph_error.push_back(massTotal+0.000001);

            f_graph.endEdit();
        }

    }


}





template <class DataTypes, class MassType>
SReal MuscleActivation<DataTypes, MassType>::getElementMass(unsigned int index) const
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    SReal mass = vertexMass[index] * massLumpingCoeff;

    return mass;
}


template <class DataTypes, class MassType>
Eigen::Vector3d MuscleActivation<DataTypes, MassType>::centerMass(const VecCoord& x, int /*numVert*/)
{	/*
			sum(m_i * r_i)		sum(r_i)
	cm = 	---------------  =	--------
				sum(m_i)		numVert
	*/

	MassVector vertexMass = vertexMassInfo.getValue();
	/*const MassVector &edgeMass = edgeMassInfo.getValue();

	double totalMass = 0.0;
	
	for(int i=0; i < _topology->getNbEdges(); i++)
	{
		unsigned idVertexA = _topology->getEdge(i)[0];
		unsigned idVertexB = _topology->getEdge(i)[1];
	
		float m = edgeMass[i] * 0.5f;

		vertexMass[idVertexA] += m;
		vertexMass[idVertexB] += m;
	}*/

	Eigen::Matrix<double, 3, 1> cm;
	cm(0) = cm(1) = cm(2) = 0.0;

	//cout << "\ngetElementMass " << getElementMass(0) << endl;
	//cout << "\nvertexMass " << vertexMass[0] << endl;
	
	for(int i=0; i < x.size(); i++)
	{
		//totalMass += vertexMass[i];

		cm(0) += x[i][0] * m_massHelper[i];
		cm(1) += x[i][1] * m_massHelper[i];
		cm(2) += x[i][2] * m_massHelper[i];
	}

	cm(0) = cm(0)*m_invTotalMass;
	cm(1) = cm(1)*m_invTotalMass;
	cm(2) = cm(2)*m_invTotalMass;

	return cm;
}

template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::multLagrange(Eigen::Vector3d &lambdaForce, Eigen::Vector3d &lambdaTorque, const VecCoord& x, Eigen::Vector3d force, Eigen::Vector3d torque, Eigen::Vector3d centerOfMass, int numVert)
{
	Eigen::Matrix3d X, X2, mat, matInv;
	Eigen::Vector3d pos;
	Eigen::Vector3d r, temp, temp2, temp3;
	temp(0) = temp(1) = temp(2) = temp2(0) = temp2(1) = temp2(2) = temp3(0) = temp3(1) = temp3(2) = 0.0;

	for(int i=0; i < x.size(); i++)
	{
		pos(0) = x[i][0];	pos(1) = x[i][1];	pos(2) = x[i][2];
		
		r = pos-centerOfMass;

		temp(0) += r(0); temp(1) += r(1); temp(2) += r(2);
		temp2(0) += r(0)*r(0); temp2(1) += r(1)*r(1); temp2(2) += r(2)*r(2);
		temp3(0) += r(0)*r(1); temp3(1) += r(1)*r(2); temp3(2) += r(0)*r(2);
	}

	X(0,0) = X(1,1) = X(2,2) = 0.0;
	X(0,1) = -temp(2);	X(0,2) = temp(1);
	X(1,0) = temp(2);	X(1,2) = -temp(0);
	X(2,0) = -temp(1);	X(2,1) = temp(0);

	// -X2
	X2(0,0) = temp2(1)+temp2(2);	X2(0,1) = -temp3(0);			X2(0,2) = -temp3(2);
	X2(1,0) = X2(0,1);				X2(1,1) = temp2(0)+temp2(2);	X2(1,2) = -temp3(1);
	X2(2,0) = X2(0,2);				X2(2,1) = X2(1,2);				X2(2,2) = temp2(0)+temp2(1);

	double invNum = double(1.0/numVert);
	mat = invNum * X * X + X2;
	//matInv = mat.inverse();
	
	//lambdaTorque = matInv * (-(1.0/numVert) * X * force + torque);
	lambdaTorque = mat.inverse() * (-invNum * X * force + torque);
	lambdaForce = invNum * (force + X*lambdaTorque);
}

template <class DataTypes, class MassType>
Eigen::Vector3d MuscleActivation<DataTypes, MassType>::createForceCorrective(Eigen::Vector3d position, Eigen::Vector3d centerOfMass, Eigen::Vector3d lambdaForce,  Eigen::Vector3d lambdaTorque)
{
	Eigen::Vector3d forceC = -lambdaForce + (position-centerOfMass).cross(lambdaTorque);

	return forceC;
}



//TODO: special case for Rigid Mass
template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::getElementMass(unsigned int index, defaulttype::BaseMatrix *m) const
{
    static const defaulttype::BaseMatrix::Index dimension = (defaulttype::BaseMatrix::Index) defaulttype::DataTypeInfo<Deriv>::size();
    if (m->rowSize() != dimension || m->colSize() != dimension) m->resize(dimension,dimension);

    m->clear();
    AddMToMatrixFunctor<Deriv,MassType>()(m, vertexMassInfo.getValue()[index] * massLumpingCoeff, 0, 1);
}

template <class DataTypes, class MassType>
void MuscleActivation<DataTypes, MassType>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    
    if (!vparams->displayFlags().getShowBehaviorModels()) return;
	if (!this->mstate) return;

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    const MassVector &vertexMass= vertexMassInfo.getValue();

    const VecCoord& x = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    Coord gravityCenter;
    Real totalMass=0.0;

	std::vector<  defaulttype::Vector3 > points;
	for (unsigned int i=0; i<x.size(); i++)
	{
		if(i < 330) continue;
		defaulttype::Vector3 p;
		p = DataTypes::getCPos(x[i]);

		points.push_back(p);
		gravityCenter += x[i]*vertexMass[i]*massLumpingCoeff;
		totalMass += vertexMass[i]*massLumpingCoeff;
	}
 
    vparams->drawTool()->drawPoints(points, 2, defaulttype::Vec<4,float>(1,0,0,1));

	for(unsigned int i = 0; i < m_orderedMuscles.size(); i++){
		for(unsigned int j = 0; j < m_orderedMuscles.at(i).size(); j++){
			for(unsigned int k = 0; k < m_orderedMuscles.at(i).at(j).size()-1; k++){
				for(unsigned int l = 0; l < m_orderedMuscles.at(i).at(j).at(k).size(); l++){
					std::vector< defaulttype::Vector3 >  points;
					defaulttype::Vector3 pa, pb;
					pa = DataTypes::getCPos(x[m_orderedMuscles[i][j][k+1][l]]);
					pb = DataTypes::getCPos(x[m_orderedMuscles[i][j][k][l]]);
					points.push_back(pa);
					points.push_back(pb);
					Vec<4,float> color2 = Vec<4,float>(0.0f , 1.0f , 0.0f,1.0f);
					vparams->drawTool()->drawLines(points,1,color2 );
				}																				  
			}
		}
	}

    if(showCenterOfGravity.getValue())
    {
        glBegin (GL_LINES);
        glColor4f (1,1,0,1);
        glPointSize(5);
        gravityCenter /= totalMass;
        for(unsigned int i=0 ; i<Coord::spatial_dimensions ; i++)
        {
            Coord v;
            v[i] = showAxisSize.getValue();
            helper::gl::glVertexT(gravityCenter-v);
            helper::gl::glVertexT(gravityCenter+v);
        }
        glEnd();
    }
#endif /* SOFA_NO_OPENGL */
}

} // namespace mass

} // namespace component

} // namespace sofa

#endif

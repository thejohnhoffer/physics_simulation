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
#define SOFA_COMPONENT_MASS_PatternPlusFileActivation_CPP
#include "PatternPlusFileActivation.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/gl/Axis.h>

namespace sofa
{

namespace component
{

namespace mass
{

using namespace sofa::defaulttype;


#ifndef SOFA_FLOAT

template <>
Vec6d PatternPlusFileActivation<Vec3dTypes, double>::getMomentum ( const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& vx, const DataVecDeriv& vv ) const
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    const MassVector &edgeMass= edgeMassInfo.getValue();

    helper::ReadAccessor< DataVecCoord > x = vx;
    helper::ReadAccessor< DataVecDeriv > v = vv;

    Vec6d momentum;
    for( unsigned int i=0 ; i<v.size() ; i++ )
    {
        Deriv linearMomentum = v[i] * vertexMass[i];
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[i], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    for( int i=0 ; i<_topology->getNbEdges() ; ++i )
    {
        unsigned v0 = _topology->getEdge(i)[0];
        unsigned v1 = _topology->getEdge(i)[1];

        // is it correct to share the edge mass between the 2 vertices?
        double m = edgeMass[i] * 0.5;

        Deriv linearMomentum = v[v0] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[v0], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];

        linearMomentum = v[v1] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        angularMomentum = cross( x[v1], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    return momentum;
}

#endif
#ifndef SOFA_DOUBLE

template <>
Vec6d PatternPlusFileActivation<Vec3fTypes, float>::getMomentum ( const core::MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& vx, const DataVecDeriv& vv ) const
{
    const MassVector &vertexMass= vertexMassInfo.getValue();
    const MassVector &edgeMass= edgeMassInfo.getValue();

    helper::ReadAccessor< DataVecCoord > x = vx;
    helper::ReadAccessor< DataVecDeriv > v = vv;

    Vec6d momentum;
    for( unsigned int i=0 ; i<v.size() ; i++ )
    {
        Deriv linearMomentum = v[i] * vertexMass[i];
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[i], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    for( int i=0 ; i<_topology->getNbEdges() ; ++i )
    {
        unsigned v0 = _topology->getEdge(i)[0];
        unsigned v1 = _topology->getEdge(i)[1];

        // is it correct to share the edge mass between the 2 vertices?
        float m = edgeMass[i] * 0.5f;

        Deriv linearMomentum = v[v0] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        Deriv angularMomentum = cross( x[v0], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];

        linearMomentum = v[v1] * m;
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[j] += linearMomentum[j];
        angularMomentum = cross( x[v1], linearMomentum );
        for( int j=0 ; j<DataTypes::spatial_dimensions ; ++j ) momentum[3+j] += angularMomentum[j];
    }

    return momentum;
}


#endif







SOFA_DECL_CLASS(PatternPlusFileActivation)

// Register in the Factory
int PatternPlusFileActivationClass = core::RegisterObject("Define a specific mass for each particle")
#ifndef SOFA_FLOAT
        .add< PatternPlusFileActivation<Vec3dTypes,double> >()
        .add< PatternPlusFileActivation<Vec2dTypes,double> >()
        .add< PatternPlusFileActivation<Vec1dTypes,double> >()
#endif
#ifndef SOFA_DOUBLE
        .add< PatternPlusFileActivation<Vec3fTypes,float> >()
        .add< PatternPlusFileActivation<Vec2fTypes,float> >()
        .add< PatternPlusFileActivation<Vec1fTypes,float> >()
#endif
        ;

/*#ifndef SOFA_FLOAT
template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<Vec3dTypes,double>;
template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<Vec2dTypes,double>;
template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<Vec1dTypes,double>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<Vec3fTypes,float>;
template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<Vec2fTypes,float>;
template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<Vec1fTypes,float>;
#endif*/

/*#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_MASS_PatternPlusFileActivation_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<defaulttype::Vec3dTypes,double>;
extern template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<defaulttype::Vec2dTypes,double>;
extern template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<defaulttype::Vec1dTypes,double>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<defaulttype::Vec3fTypes,float>;
extern template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<defaulttype::Vec2fTypes,float>;
extern template class SOFA_SiElegansPlugin_API PatternPlusFileActivation<defaulttype::Vec1fTypes,float>;
#endif
#endif*/


} // namespace mass

} // namespace component

} // namespace sofa


#define SOFA_COMPONENT_FORCEFIELD_Friction_CPP

#include "Friction.inl"

#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ForceField.inl>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(Friction)

// Register in the Factory
int FrictionClass = core::RegisterObject("Generic Tetrahedral finite elements")
#ifndef SOFA_FLOAT
.add< Friction<Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< Friction<Vec3fTypes> >()
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_SiElegansPlugin_API Friction<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_SiElegansPlugin_API Friction<Vec3fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

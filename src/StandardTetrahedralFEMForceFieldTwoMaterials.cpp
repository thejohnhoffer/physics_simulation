#define SOFA_COMPONENT_FORCEFIELD_MYSTANDARDTETRAHEDRALFEMFORCEFIELD_CPP

#include "StandardTetrahedralFEMForceFieldTwoMaterials.inl"

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

SOFA_DECL_CLASS(StandardTetrahedralFEMForceFieldTwoMaterials)

// Register in the Factory
int MyStandardTetrahedralFEMForceFieldClass = core::RegisterObject("Generic Tetrahedral finite elements")
#ifndef SOFA_FLOAT
.add< StandardTetrahedralFEMForceFieldTwoMaterials<Vec3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< StandardTetrahedralFEMForceFieldTwoMaterials<Vec3fTypes> >()
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_SiElegansPlugin_API StandardTetrahedralFEMForceFieldTwoMaterials<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_SiElegansPlugin_API StandardTetrahedralFEMForceFieldTwoMaterials<Vec3fTypes>;
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

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
#include "StimuliBarycentricPenalityContact.inl"
#include <sofa/core/ObjectFactory.h>
#include <SofaMeshCollision/BarycentricContactMapper.h>
#include <SofaMeshCollision/RigidContactMapper.h>
#include <SofaBaseCollision/RigidCapsuleModel.h>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace core::collision;
using simulation::Node;

SOFA_DECL_CLASS(StimuliBarycentricPenalityContact)

Creator<Contact::Factory, StimuliBarycentricPenalityContact<SphereModel, SphereModel> > SphereSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<SphereModel, RigidSphereModel> > SphereRigidSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidSphereModel, RigidSphereModel> > RigidSphereRigidSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<SphereModel, PointModel> > SpherePointPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidSphereModel, PointModel> > RigidSpherePointPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<PointModel, PointModel> > PointPointPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<LineModel, PointModel> > LinePointPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<LineModel, LineModel> > LineLinePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<LineModel, SphereModel> > LineSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<LineModel, RigidSphereModel> > LineRigidSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<TriangleModel, SphereModel> > TriangleSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<TriangleModel, RigidSphereModel> > TriangleRigidSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<TriangleModel, PointModel> > TrianglePointPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<TriangleModel, LineModel> > TriangleLinePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<TriangleModel, TriangleModel> > TriangleTrianglePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, TriangleModel> > CapsuleTrianglePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, LineModel> > CapsuleLinePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, CapsuleModel> > CapsuleCapsulePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, SphereModel> > CapsuleSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, RigidSphereModel> > CapsuleRigidSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<OBBModel, OBBModel> > OBBOBBPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, OBBModel> > CapsuleOBBPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<SphereModel, OBBModel> > SphereOBBPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidSphereModel, OBBModel> > RigidSphereOBBPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<TriangleModel, OBBModel> > TriangleOBBPenalityContactClass("SiElegans",true);

Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidCapsuleModel, TriangleModel> > RigidCapsuleTrianglePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidCapsuleModel, LineModel> > RigidCapsuleLinePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidCapsuleModel, RigidCapsuleModel> > RigidCapsuleRigidCapsulePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, RigidCapsuleModel> > CapsuleRigidCapsulePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidCapsuleModel, SphereModel> > RigidCapsuleSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidCapsuleModel, RigidSphereModel> > RigidCapsuleRigidSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<RigidCapsuleModel, OBBModel> > RigidCapsuleOBBPenalityContactClass("SiElegans",true);

Creator<Contact::Factory, StimuliBarycentricPenalityContact<CylinderModel, CylinderModel> > CylinderCylinderPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CylinderModel, TriangleModel> > CylinderTrianglePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CylinderModel, RigidCapsuleModel> > CylinderRigidCapsulePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CapsuleModel, CylinderModel> > CapsuleCylinderPenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CylinderModel, SphereModel> > CylinderSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CylinderModel, RigidSphereModel> > CylinderRigidSpherePenalityContactClass("SiElegans",true);
Creator<Contact::Factory, StimuliBarycentricPenalityContact<CylinderModel, OBBModel> > CylinderOBBPenalityContactClass("SiElegans",true);
} // namespace collision

} // namespace component

} // namespace sofa



#ifndef SOFA_COMPONENT_FORCEFIELD_SiElegansAnimationLoop_H
#define SOFA_COMPONENT_FORCEFIELD_SiElegansAnimationLoop_H

#if !defined(__GNUC__) || (__GNUC__ > 3 || (_GNUC__ == 3 && __GNUC_MINOR__ > 3))
#pragma once
#endif

#include "initPlugin.h"
#include "Stimuli.h"
#include "SaveScreenShots.h"
#include "InterfaceManagerActivation.h"

#include <sofa/simulation/common/DefaultAnimationLoop.h>


namespace sofa
{

namespace simulation
{



//***************** Tetrahedron FEM code for several elastic models: StandardTetrahedralFEMForceField*******************************************************************
//********************************** Based on classical discretization : Fi=-Bi^T S V and Kij=Bi^T N Bj +Di^T S Dj **********************************************
//***************************************** where Bi is the strain displacement (6*3 matrix), S SPK tensor N=dS/dC, Di shape vector ************************************
//**************************** Code dependant on HyperelasticMatrialFEM and inherited classes *********************************************************************

/** Compute Finite Element forces based on tetrahedral elements.
*/
class SOFA_SiElegansPlugin_API SiElegansAnimationLoop: public sofa::simulation::DefaultAnimationLoop
{
  public:

    typedef sofa::simulation::DefaultAnimationLoop Inherit;
    typedef sofa::core::objectmodel::BaseContext BaseContext;
    typedef sofa::core::objectmodel::BaseObjectDescription BaseObjectDescription;
    SOFA_CLASS(SiElegansAnimationLoop,sofa::simulation::DefaultAnimationLoop);
protected:
    SiElegansAnimationLoop(simulation::Node* gnode = NULL);

    virtual ~SiElegansAnimationLoop();
public:

    /// Set the simulation node to the local context if not specified previously
    virtual void init();

    /// perform one animation step
    virtual void step(const core::ExecParams* params, SReal dt);

	Data<std::string> f_getURL;
	Data<std::string> f_postURL;
	Data<std::string> f_finishURL;

	Data<int> f_finishingTime;
	Data<int> f_syncTime;
	float m_currentSyncTime;

	std::vector<component::behaviormodel::Stimuli*> m_stimuliVec;
	sofa::component::misc::SaveScreenShots* m_screensaver;
	sofa::component::mass::InterfaceManagerActivation<defaulttype::Vec3dTypes,double>* m_muscleActivation;

	Data<bool> f_writeTraces;

	std::vector<std::string> &split(std::string &s, char delim, std::vector<std::string> &elems);

	std::vector<std::string> split(std::string &s, char delim);

};

} // namespace component

} // namespace sofa

#endif 

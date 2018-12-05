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
*                    Si Elegans :: Behavioural models  :: [Stimulus]           *
*                                                                             *
* Contact information: info@vicotmech.org                                     *
******************************************************************************/
#include "[Stimulus].h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/helper/gl/template.h>
#include <sofa/core/ObjectFactory.h>
#include <iostream>
#include <string.h>
#include <sofa/defaulttype/BoundingBox.h>
#include <SofaBaseMechanics/MechanicalObject.h>

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/rapidjson.h"
#include <iostream>
using namespace rapidjson;

namespace sofa
{

namespace component
{

namespace behaviormodel
{


SOFA_DECL_CLASS([Stimulus])

int [Stimulus]Class = core::RegisterObject("[Stimulus] class")
        .add< [Stimulus] >()
        .addLicense("LGPL")
        .addAuthor("PL")
        ;

[Stimulus]::[Stimulus](): Stimuli()
{
}

[Stimulus]::~[Stimulus]()
{
   
}

void [Stimulus]::init()
{
	DOUT("init [Stimulus] -> "  << this->getName());
	Stimuli::init();

}

void [Stimulus]::reset() // resets the model
{
      resetStimuli();
      // TODO init stimuli for dynamic ones
}

void [Stimulus]::step(){
	//sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
	//sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());

	DOUT("step() - [Stimulus] -> " << this->getName());
        // TODO update the stimuli and generate stimuli values -> the main part of the stimulus-ception
}

void [Stimulus]::updatePosition(double /*dt*/)
{
}

void [Stimulus]::draw(const core::visual::VisualParams* /*vparams*/)
{
}

void [Stimulus]::exportOBJ(std::string /*name*/, std::ostream* /*out*/, std::ostream* /*mtl*/, int& /*vindex*/, int& /*nindex*/, int& /*tindex*/)
{
}

void [Stimulus]::updateVisual()
{
}

void [Stimulus]::computeBBox(const core::ExecParams*  /*params*/ )
{
}

void [Stimulus]::createStimuliString(){
	Stimuli::createStimuliString();
	DOUT("createStimuliString() - [Stimulus] -> " << this->getName());
	//include first timestep
	std::ostringstream out;
	
	// out << ... TODO  format the stimuli values to the required output string;
	// out << "#" << m_neuronIDmap["ADAL"] << "," << 65522 << "," << 0.564;			
	
	this->m_stimuliString = out.str();
  DOUT("--- output stimuli str:" << this->m_stimuliString);
}

void [Stimulus]::resetStimuli(){
	DOUT("resetStimuli() - [Stimulus] -> " << this->getName());
  
    // TODO reset the internal stimuli values
}

} // namespace behaviormodel

} // namespace component

} // namespace sofa


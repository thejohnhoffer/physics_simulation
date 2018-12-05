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
*                    Si Elegans :: Behavioural models  :: StimuliMechano           *
*                                                                             *
* Contact information: info@vicotmech.org                                     *
******************************************************************************/
#include "StimuliMechano.h"
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


SOFA_DECL_CLASS(StimuliMechano)

int StimuliMechanoClass = core::RegisterObject("StimuliMechano class")
        .add< StimuliMechano >()
        .addLicense("LGPL")
        .addAuthor("PL")
        ;

StimuliMechano::StimuliMechano(): Stimuli(),
   f_fileNodeNeuronDist(initData(&f_fileNodeNeuronDist,"file_NodeNeuronDist","file with distances between nodes and neurons"))
{
   f_fileNodeNeuronDist.setGroup("config files");
}

StimuliMechano::~StimuliMechano()
{
   
}

void StimuliMechano::init()
{
  DOUT("init StimuliMechano "  << this->getName());
  Stimuli::init();
  
  COUT("loadNodeNeuronDistances");
  //load the file that contains the distances between nodes and neurons
  if (!f_fileNodeNeuronDist.getValue().empty())
  {
    m_nodeNeuronDistances = this->loadNodeNeuronDistances(f_fileNodeNeuronDist.getFullPath().c_str());
  }
  
  //m_mechanoForces.resize(m_neuronIDmap.size(), 0.0);
}

void StimuliMechano::reset() // resets the model
{
      resetStimuli();
}

void StimuliMechano::step(){
  //sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  //sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());

  // TODO update the stimuli and generate stimuli values -> the main part of the stimulus-ception
  // -> currently the force stimuli are added by other components via the addMechanoStimuli method, upon contact detection
}

void StimuliMechano::updatePosition(double /*dt*/)
{
}

void StimuliMechano::draw(const core::visual::VisualParams* /*vparams*/)
{
}

void StimuliMechano::exportOBJ(std::string /*name*/, std::ostream* /*out*/, std::ostream* /*mtl*/, int& /*vindex*/, int& /*nindex*/, int& /*tindex*/)
{
}

void StimuliMechano::updateVisual()
{
}

void StimuliMechano::computeBBox(const core::ExecParams*  /*params*/ )
{
}

void StimuliMechano::createStimuliString(){
  Stimuli::createStimuliString();
  //include first timestep
  std::ostringstream strs;
  
  DOUT("createStimuliString() - mechano " << this->getName());
  for(t_mechanoForces::iterator it=m_mechanoForces.begin(); it!=m_mechanoForces.end(); ++it)
      strs << "#" << this->m_neuronIDmap[it->first] << ",65524," << it->second;
 
  this->m_stimuliString = strs.str();
  DOUT("--- output stimuli str:" << this->m_stimuliString);

  std::cout << "####" << std::endl;
  std::cout << "BDUL " << m_mechanoForces["BDUL"] << std::endl;
  std::cout << "BDUR " << m_mechanoForces["BDUR"] << std::endl;
  std::cout << "LUAL " << m_mechanoForces["LUAL"] << std::endl;
  std::cout << "LUAR " << m_mechanoForces["LUAR"] << std::endl;
}

void StimuliMechano::resetStimuli(){
  
  DOUT("resetStimuli() - mechano " << this->getName()); 
  if (m_mechanoForces.empty()) return;
  t_mechanoForces::iterator pit=m_mechanoForces.begin();
  t_mechanoForces::iterator it=m_mechanoForces.begin();++it;
  for(; it!=m_mechanoForces.end(); ++it){
      if (pit->second <=0.0)
       m_mechanoForces.erase(pit);
      else
       pit->second = 0.0;
      pit = it;
  }
      if (pit->second <=0.0)
       m_mechanoForces.erase(pit);
      else
       pit->second = 0.0;
  
  //m_mechanoForces.clear();
}


void StimuliMechano::addMechanoStimuli(unsigned int nodeIndex, double forceValue){
  DOUT("addMechanoStimuli() - mechano " << this->getName());
  // follow function get_cuticle_node_to_neuron_list_and_weight(nearestNode, touchType) in predefined_stimuli.py
  
  float g_weight = 1.0; // gaussian weight ... TODO adapt to consider the force smoothing by gaussians -> currently no contact smoothing considered
  double worm_length = this->m_minMaxXYZ[1]-this->m_minMaxXYZ[0];
  double worm_width = this->m_minMaxXYZ[5]-this->m_minMaxXYZ[4];
 
  for(t_neuronDistances::iterator it = m_nodeNeuronDistances[nodeIndex].begin(); it != m_nodeNeuronDistances[nodeIndex].end(); ++it)
  {
        // python code (predefined_stimuli.py) follows
        // if (dist_to_nodes[node] < 2.0 * WORM_WIDTH / WORM_LENGHT): #FIXME play with limitation on ditance, e.g. 3*sigma of the gaussian ? Currently 2 x worm width
        //   # weight = weight + 1.0/(1+dist_to_nodes[node]*10) * g_weight # gaussian * 1/dist considered ... FIXME change here for other stress 'diffusion' law
        //   weight = max(weight,1.0/(1+dist_to_nodes[node]*10) * g_weight) # FIXME use different force accumulation law  
	if (it->second < 2.0 * worm_width/worm_length){// limit the maximum distance to double of the worm's width (remember distances are normalized to 1=WORM_LENGTH)
          //m_mechanoForces.at[it->first] += it->second*forceValue;
	  double weight = 1.0/(it->second*1000) * g_weight;
	  if (weight > 0.0001) // FIXME only to limit the data flow
          {
             t_mechanoForces::iterator itf = m_mechanoForces.find(it->first);
             if (itf != m_mechanoForces.end())
                itf->second += weight * forceValue;
	     else 
                m_mechanoForces[it->first] = weight * forceValue;
          }
        }
  }
}



} // namespace behaviormodel

} // namespace component

} // namespace sofa


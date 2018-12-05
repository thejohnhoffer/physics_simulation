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
*                    Si Elegans :: Behavioural models  :: StimuliThermo           *
*                                                                             *
* Contact information: info@vicotmech.org                                     *
******************************************************************************/
#include "StimuliThermo.h"
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
#include <cmath>
using namespace rapidjson;

namespace sofa
{

namespace component
{

namespace behaviormodel
{


SOFA_DECL_CLASS(StimuliThermo)

int StimuliThermoClass = core::RegisterObject("StimuliThermo class")
        .add< StimuliThermo >()
        .addLicense("LGPL")
        .addAuthor("PL")
        ;

StimuliThermo::StimuliThermo(): Stimuli()
   , f_fileNeuronSomaPos(initData(&f_fileNeuronSomaPos,"file_NeuronSomaPos","file defining the Neuron soma positions on the A-P axis"))
   , f_envTemperature(initData(&f_envTemperature, (float)21.0, "envTemperature", "environment temperature"))
   , f_heatDiffusion(initData(&f_heatDiffusion, (float)1.0, "heatDiffusion", "heat diffusion constant"))
   , f_gradientOn(initData(&f_gradientOn, false, "gradientOn", "turn on Left to Right temperature gradient"))
   , f_LRgradDistance(initData(&f_LRgradDistance, (float)200.0, "LRgradDistance", "distance between the left and right extrema of the temperature gradient"))
   , f_gradLeft(initData(&f_gradLeft, (float)21.0, "heatLeft", "heat value on the left side of the dish"))
   , f_gradRight(initData(&f_gradRight, (float)21.0, "heatRight", "heat value on the right side of the dish"))
   , f_heatPoints(initData(&f_heatPoints, "heatPoints", "list of static temperature heat points"))
   , f_heatPointTemperature(initData(&f_heatPointTemperature, "heatPointTemperature", "list of temperatures for the heat points"))
   , f_heatPointTiming(initData(&f_heatPointTiming, "heatPointTiming", "list of heat points activation time"))
   , f_heatPointDuration(initData(&f_heatPointDuration, "heatPointDuration", "list of heat points active duration"))
{
  f_fileNeuronSomaPos.setGroup("config files");

  f_envTemperature.setGroup("Thermo Environment");
  f_heatDiffusion.setGroup("Thermo Environment");

  f_LRgradDistance.setGroup("Thermo L-R Gradient");
  f_gradLeft.setGroup("Thermo L-R Gradient");
  f_gradRight.setGroup("Thermo L-R Gradient");

  f_heatPoints.setGroup("Thermo Heat Point");
  f_heatPointTemperature.setGroup("Thermo Heat Point");
  f_heatPointTiming.setGroup("Thermo Heat Point");
  f_heatPointDuration.setGroup("Thermo Heat Point");
}

StimuliThermo::~StimuliThermo()
{
   
}

void StimuliThermo::init()
{
  DOUT("init StimuliThermo " << this->getName());
  Stimuli::init();

  COUT("loadNeuronSomaPositions");
  //load the file that contains the positions of the neuron's soma on the Anterior-to-posterios axis (normalized to 1)
  if (!f_fileNeuronSomaPos.getValue().empty())
  {
    m_neuronSomaPos = this->loadNeuronPositions(f_fileNeuronSomaPos.getFullPath().c_str());
  }
  
  m_envTemperature = f_envTemperature.getValue();
  m_heatDiffusion = f_heatDiffusion.getValue();
  
  m_heatGradientOn = f_gradientOn.getValue();
  if (m_heatGradientOn){
    m_LRgradDistance = f_LRgradDistance.getValue();
    m_gradLeft = f_gradLeft.getValue();
    m_gradRight = f_gradRight.getValue();
  }
  
  m_heatPoints = f_heatPoints.getValue();
  m_heatPointTemperature = f_heatPointTemperature.getValue();
  m_heatPointTiming = f_heatPointTiming.getValue();
  m_heatPointDuration = f_heatPointDuration.getValue();
}

void StimuliThermo::reset() // resets the model
{
      resetStimuli();
      // TODO init stimuli for dynamic ones
}

float StimuliThermo::getTemepratureStatic2(double r2, float t0, float sigma)// distance^2, steady temperature in the middle, heat_diffusion i.e. sigma
{
    // static: non-normalized gaussian used t = t0 * exp(- (r^2)/(2*sigma^2))    
    double tmp = r2/(sigma*sigma);
    return t0 * std::exp(- 0.5*tmp);
}

void StimuliThermo::step(){
  //sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  //sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());

  
  DOUT("StimuliThermo compute stimuli " << this->getName());
  for (t_neuronPos::iterator it=m_neuronSomaPos.begin(); it != m_neuronSomaPos.end(); ++it)
  {
      std::string neuronName = it->first;
      float _pos = it->second;
      defaulttype::Vec3dTypes::Coord pos = WormRelPosToWorldPos(_pos, 0.5, 0.5);// consider them in the middle // TODO discern left/right
      //DOUT( it->first << "-> x: " << pos[0] << " y: " << pos[1] << " z: " << pos[2] );
      
      // gradient left - right
      if (m_heatGradientOn){
          //DOUT("StimuliThermo  compute Left to Right gradient values");
          if (m_writeTraces) std::cout << "g";

          float temp = m_gradLeft + ( m_gradRight - m_gradLeft)/m_LRgradDistance * (pos[0] + m_LRgradDistance/2.0);
          if (pos[0] < -m_LRgradDistance/2.0) temp = m_gradLeft;
          if (pos[0] > m_LRgradDistance/2.0) temp = m_gradRight;
	  DOUT ("relDist " << (pos[0] + m_LRgradDistance/2.0)/m_LRgradDistance << " temp " << temp)
		

          t_temepratureSensoryOut::iterator outIt= m_neuronTemepratureSensoryOut.find(neuronName); 
          if (outIt == m_neuronTemepratureSensoryOut.end())
            m_neuronTemepratureSensoryOut[neuronName] = temp; // TODO
          else           
            m_neuronTemepratureSensoryOut[neuronName] = temp; // TODO how to combine temperatures form different sources ? use max ? average ? 
      }
      
      // heat points 
      if (!m_heatPoints.empty()){
        //DOUT("StimuliThermo  compute heat point values");
        if (m_writeTraces) std::cout << ".";
      }
      for (unsigned int i=0; i< m_heatPoints.size(); i++)
      {
        float relTime = this->getContext()->getTime() - m_heatPointTiming[i]/1000;        // sofa timing in sec, our in msec
        if ((relTime >= 0) && (relTime <= m_heatPointDuration[i]/1000))
        {
            double distX = pos[0] - m_heatPoints[i][0];
            double distZ = pos[2] - m_heatPoints[i][1];
            double dist = distX*distX + distZ*distZ;
	    double temp = getTemepratureStatic2(dist, m_heatPointTemperature[i], m_heatDiffusion);
	    DOUT ("dist " << dist << " temp " << temp)
            
            t_temepratureSensoryOut::iterator outIt= m_neuronTemepratureSensoryOut.find(neuronName);  
            if (outIt == m_neuronTemepratureSensoryOut.end())
              m_neuronTemepratureSensoryOut[neuronName] = temp; // TODO
            else           
              m_neuronTemepratureSensoryOut[neuronName] = temp; // TODO how to combine temperatures form different sources ? use max ? average ? 
        }      
      }        
      
      //DOUT("StimuliThermo  special neuron handling");
      //if (m_writeTraces) std::cout << ".";
      // TODO special handling for AFD (amphid), AWC (amphid), FLP (branched, hot shock), PVD (branched, cold shock), and PHC (tail)
      
  }
  
}

void StimuliThermo::updatePosition(double /*dt*/)
{
}

void StimuliThermo::draw(const core::visual::VisualParams* /*vparams*/)
{
}

void StimuliThermo::exportOBJ(std::string /*name*/, std::ostream* /*out*/, std::ostream* /*mtl*/, int& /*vindex*/, int& /*nindex*/, int& /*tindex*/)
{
}

void StimuliThermo::updateVisual()
{
}

void StimuliThermo::computeBBox(const core::ExecParams*  /*params*/ )
{
}

void StimuliThermo::createStimuliString(){
  Stimuli::createStimuliString();
  DOUT("createStimuliString() - StimuliThermo " << this->getName());
  //include first timestep
  std::ostringstream out;
  
  // out << ... TODO  format the stimuli values to the required output string;
  // out << "#" << m_neuronIDmap["ADAL"] << "," << 65522 << "," << 0.564;     
  for (std::map< std::string, float >::iterator it = m_neuronTemepratureSensoryOut.begin(); it != m_neuronTemepratureSensoryOut.end(); ++it)
  { // neuron temperature sensing
      out << "#" << m_neuronIDmap[it->first] << "," << 65525 << "," << it->second;      
  }

  this->m_stimuliString = out.str();
  DOUT("--- output stimuli str:" << this->m_stimuliString);
}

void StimuliThermo::resetStimuli(){
  DOUT("resetStimuli() - StimuliThermo " << this->getName());
  
  for (std::map< std::string, float >::iterator it = m_neuronTemepratureSensoryOut.begin(); it != m_neuronTemepratureSensoryOut.end(); ++it)
    it->second = m_envTemperature;
}

} // namespace behaviormodel

} // namespace component

} // namespace sofa


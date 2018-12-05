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
*                    Si Elegans :: Behavioural models  :: StimuliProprio           *
*                                                                             *
* Contact information: info@vicotmech.org                                     *
******************************************************************************/
#include "StimuliProprio.h"
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

//#define PI   3.14159265358979323846

using namespace rapidjson;

namespace sofa
{

namespace component
{

namespace behaviormodel
{


SOFA_DECL_CLASS(StimuliProprio)

int StimuliProprioClass = core::RegisterObject("StimuliProprio class")
        .add< StimuliProprio >()
        .addLicense("LGPL")
        .addAuthor("PL")
        ;

StimuliProprio::StimuliProprio(): Stimuli(),
   f_fileNeuronSomaPos(initData(&f_fileNeuronSomaPos,"file_NeuronSomaPos","file defining the Neuron soma positions on the A-P axis")),
   f_fileMuscleExtremes(initData(&f_fileMuscleExtremes,"file_MuscleExtremes","file with node indices of muscle extremes")),
   f_fileMusclePositions(initData(&f_fileMusclePositions,"file_MusclePositions","file defining the Muscle positions wrt the 3D model (8 octants of 12 elements)")),
   f_fileMuscleSensoryConf(initData(&f_fileMuscleSensoryConf,"file_MuscleSensoryConf","file defining the Muscle to Neuron connections for implementing muscle sensation"))
{
  f_fileNeuronSomaPos.setGroup("config files");
  f_fileMuscleSensoryConf.setGroup("config files");
  f_fileMusclePositions.setGroup("config files");
  f_fileMuscleExtremes.setGroup("config files");
}

StimuliProprio::~StimuliProprio()
{
   
}

void StimuliProprio::init()
{
  DOUT("init StimuliProprio " << this->getName());
  Stimuli::init();

  COUT("loadNeuronSomaPositions");
  //load the file that contains the positions of the neuron's soma on the Anterior-to-posterios axis (normalized to 1)
  if (!f_fileNeuronSomaPos.getValue().empty())
  {
    m_neuronSomaPos = this->loadNeuronPositions(f_fileNeuronSomaPos.getFullPath().c_str());
  }
    
  //load the file that contains the indices of the nodes of the extremes of muscles
  COUT("loadMuscleExtremes");
  if (!f_fileMuscleExtremes.getValue().empty())
  {
    m_extremeIndices = this->loadMuscleExtremes(f_fileMuscleExtremes.getFullPath().c_str());
  }

  COUT("loadMusclePosition");
  //load the file that contains the configuration of the muscle to neuron connections for muscle based proprioception
  if (!f_fileMusclePositions.getValue().empty())
  {
    m_musclePositionConf = this->loadMuscleConfiguration(f_fileMusclePositions.getFullPath().c_str());
  }
  
  COUT("loadMuscleSensoryConfiguration");
  //load the file that contains the configuration of the muscle to neuron connections for muscle based proprioception
  if (!f_fileMuscleSensoryConf.getValue().empty())
  {
    m_muscleSensoryConf = this->loadMuscleSensoryConfiguration(f_fileMuscleSensoryConf.getFullPath().c_str());
  }


  for (t_muscleSensoryConf::iterator it = m_muscleSensoryConf.begin(); it != m_muscleSensoryConf.end(); ++it)
      m_neuronMuscleSensoryOut[it->first] = 0.0;

  DOUT("computeMuscleEnds");
  m_restLengths = this->computeMuscleEnds( m_extremeIndices );
}

void StimuliProprio::reset() // resets the model
{
      resetStimuli();
      // TODO init stimuli for dynamic ones
}

// return the relative length of the muscle, wrt. restlength
float StimuliProprio::calculateMuscleRelLenght(std::string muscleName, const t_muscleEnds& initEnds, const t_muscleEnds& currentEnds)
{
    t_muscleConf mc = m_musclePositionConf[muscleName];
    float relLen =  ((currentEnds[mc.octant][mc.element][0]-currentEnds[mc.octant][mc.element][1])).norm();
    relLen /= (initEnds[mc.octant][mc.element][0]-initEnds[mc.octant][mc.element][1]).norm();
    return relLen;
}

float StimuliProprio::calculateFlectionAt( std::string neuronName, t_pointVec centralPoints)
{
  // .. get nearest of the central points and calculate angle
    float neuron_pos = m_neuronSomaPos[neuronName];
    int index = (int) floor( neuron_pos / (1.0/centralPoints.size()) + 0.5);
    if ((index >= (int) centralPoints.size() -1) || (index <= 0)) 
        return PI; // extremal points, no flection (no info)

    sofa::defaulttype::Vec< 3, float > p1,p2,p3;
    p1 = centralPoints[index-1];
    p2 = centralPoints[index];
    p3 = centralPoints[index+1];
    
    sofa::defaulttype::Vec< 3, float > v1 = p1 - p2;
    sofa::defaulttype::Vec< 3, float > v2 = p3 - p2;
	v1[1] = 0;
	v2[1] = 0;
    v1.normalize();
    v2.normalize();
	float dir = v1.cross(v2)[1];
    
    float dot = v1*v2;
	float angle;
	if(dir >= 0){
		angle = std::acos(dot);
	}
	else{
		angle = 2*PI - std::acos(dot);
	}

	//std::cout << "dir: " << dir << " " << angle << "\n";
      
    return angle;
}

void StimuliProprio::step(){
  DOUT("compute stimuli - StimuliProprio " << this->getName());

  //sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  //sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());

  // compute muscle stretch
  t_muscleEnds currentEnds = computeMuscleEnds( m_extremeIndices );

  for (t_muscleSensoryConf::iterator it = m_muscleSensoryConf.begin(); it != m_muscleSensoryConf.end(); ++it)
  {
    for (unsigned int i=0; i< it->second.size(); i++){
        // TODO += is used due to multiple inputs; simple summatory rule for sensory composition is considered
        m_neuronMuscleSensoryOut[it->first] += calculateMuscleRelLenght(it->second[i].name, m_restLengths, currentEnds) * it->second[i].weight; 

		//if(it->first == "AS1" && i == 0) std::cout << it->first << "a: " << calculateMuscleRelLenght(it->second[i].name, m_restLengths, currentEnds) << std::endl;
		//if(it->first == "VD2" && i == 1) std::cout << it->first << ": " << calculateMuscleRelLenght(it->second[i].name, m_restLengths, currentEnds) << std::endl;
	} 
  }

  /*std::cout <<  "AS1t: " << m_neuronMuscleSensoryOut["AS1"] << std::endl;
  std::cout <<  "DD2: " << m_neuronMuscleSensoryOut["DD2"] << std::endl;
  std::cout <<  "VD4: " << m_neuronMuscleSensoryOut["VD4"] << std::endl;*/
  // compute body bend
  t_pointVec centralPoints = computeCentralPoints();
  for (std::map<std::string, int>::iterator it=m_neuronIDmap.begin(); it != m_neuronIDmap.end(); ++it){
	  m_neuronBendSensoryOut[it->first] = calculateFlectionAt( it->first, centralPoints);
	  if(it->first == "VA6" || it->first == "AS3"){
		  std::cout << it->first << " " << m_neuronBendSensoryOut[it->first] << std::endl;
	  }
  }
}

void StimuliProprio::updatePosition(double /*dt*/)
{
}

void StimuliProprio::draw(const core::visual::VisualParams* /*vparams*/)
{
}

void StimuliProprio::exportOBJ(std::string /*name*/, std::ostream* /*out*/, std::ostream* /*mtl*/, int& /*vindex*/, int& /*nindex*/, int& /*tindex*/)
{
}

void StimuliProprio::updateVisual()
{
}

void StimuliProprio::computeBBox(const core::ExecParams*  /*params*/ )
{
}

void StimuliProprio::createStimuliString(){
  Stimuli::createStimuliString();
  DOUT("createStimuliString() - StimuliProprio " << this->getName());
  //include first timestep
  std::ostringstream out;
  
  for (std::map< std::string, float >::iterator it = m_neuronMuscleSensoryOut.begin(); it != m_neuronMuscleSensoryOut.end(); ++it)
  { // muscle stretch output
    out << "#" << m_neuronIDmap[it->first] << "," << 65522 << "," << it->second;      
  }
  
  for (std::map< std::string, float >::iterator it = m_neuronBendSensoryOut.begin(); it != m_neuronBendSensoryOut.end(); ++it)
  { // body bend sensing
    out << "#" << m_neuronIDmap[it->first] << "," << 65521 << "," << it->second;      
  }

  this->m_stimuliString = out.str();
  DOUT("--- output stimuli str:" << this->m_stimuliString);
}

void StimuliProprio::resetStimuli(){
  DOUT("resetStimuli() - StimuliProprio " << this->getName());
  
  for (std::map< std::string, float >::iterator it = m_neuronMuscleSensoryOut.begin(); it != m_neuronMuscleSensoryOut.end(); ++it)
    it->second = 0.0;
  for (std::map< std::string, float >::iterator it = m_neuronBendSensoryOut.begin(); it != m_neuronBendSensoryOut.end(); ++it)
    it->second = 0.0;
}

} // namespace behaviormodel

} // namespace component

} // namespace sofa


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
*                    Si Elegans :: Behavioural models  :: StimuliChemo           *
*                                                                             *
* Contact information: info@vicotmech.org                                     *
******************************************************************************/
#include "StimuliChemo.h"
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

using namespace sofa::defaulttype;

namespace sofa
{

namespace component
{

namespace behaviormodel
{


SOFA_DECL_CLASS(StimuliChemo)
float StimuliChemo::t_chemoOutputValues::dummy = -1.0;

int StimuliChemoClass = core::RegisterObject("StimuliChemo class")
        .add< StimuliChemo >()
        .addLicense("LGPL")
        .addAuthor("PL")
        ;

StimuliChemo::StimuliChemo(): Stimuli()
  , f_fileSensoryPos(initData(&f_fileSensoryPos,"file_SensoryPos","file defining the Sensory Positions"))
  , f_fileChemicalProp(initData(&f_fileChemicalProp,"file_ChemicalProp","file defining the Chemical Properties (diffusivity, time_diffusion, time_evaporation)"))
    
  , f_barrierWidth ( initData(&f_barrierWidth, (float)2.0, "barrierWidth", "barrier width relative to WORM size") )
  , f_barrierChemical(initData(&f_barrierChemical, "barrierChemical", "barrier chemical"))
  , f_barrierConcentration ( initData(&f_barrierConcentration, (float)0.0, "barrierConcentration", "barrier concentration") )
  , f_barrierVerticalOn( initData(&f_barrierVerticalOn, false, "barrierVerticalOn","sets vertical barrier on/off") )
  , f_barrierHorizontalOn( initData(&f_barrierHorizontalOn, false, "barrierHorizontalOn","sets horizontal barrier on/off") )
  , f_barrierRingOn( initData(&f_barrierRingOn, false, "barrierRingOn","sets osmotic ring barrier on/off") )
  , f_ringInnerRadius( initData(&f_ringInnerRadius, (float)10.0, "ringInnerRadius","sets osmotic ring barrier inner radius") )
    
  , f_q1Chemical(initData(&f_q1Chemical, "q1Chemical", "q1 Chemical"))
  , f_q2Chemical(initData(&f_q2Chemical, "q2Chemical", "q2 Chemical"))
  , f_q3Chemical(initData(&f_q3Chemical, "q3Chemical", "q3 Chemical"))
  , f_q4Chemical(initData(&f_q4Chemical, "q4Chemical", "q4 Chemical"))
  , f_q1Concentration ( initData(&f_q1Concentration, (float)0.0, "q1Concentration", "q1 Concentration") )
  , f_q2Concentration ( initData(&f_q2Concentration, (float)0.0, "q2Concentration", "q2 Concentration") )
  , f_q3Concentration ( initData(&f_q3Concentration, (float)0.0, "q3Concentration", "q3 Concentration") )
  , f_q4Concentration ( initData(&f_q4Concentration, (float)0.0, "q4Concentration", "q4 Concentration") )

  , f_staticPoints(initData(&f_staticPoints, "staticPoints", "static points list"))
  , f_staticChemicals(initData(&f_staticChemicals, "staticChemicals", "static chemical list"))
  , f_staticConcentrations(initData(&f_staticConcentrations, "staticConcentrations", "static concentrations list"))
  
  , f_dynamicPoints(initData(&f_dynamicPoints, "dynamicPoints", "dynamic Points list"))
  , f_dynamicChemicals(initData(&f_dynamicChemicals, "dynamicChemicals", "dynamic Chemical list"))
  , f_dynamicConcentrations(initData(&f_dynamicConcentrations, "dynamicConcentrations", "dynamic Concentrations list"))
  , f_dynamicTiming(initData(&f_dynamicTiming, "dynamicTiming", "dynamic Timing"))  

  , f_dropChemicals(initData(&f_dropChemicals, "dropChemicals", "drop Chemical list"))
  , f_dropVolumes(initData(&f_dropVolumes, "dropVolumes", "drop volume list"))
  , f_dropConcentrations(initData(&f_dropConcentrations, "dropConcentrations", "drop Concentrations list"))
  , f_dropTiming(initData(&f_dropTiming, "dropTiming", "drop Timing"))  
{
  f_fileSensoryPos.setGroup("config files");
  f_fileChemicalProp.setGroup("config files");

  f_barrierWidth.setGroup("Chemo Barrier");
  f_barrierChemical.setGroup("Chemo Barrier");
  f_barrierConcentration.setGroup("Chemo Barrier");
  f_barrierVerticalOn.setGroup("Chemo Barrier");
  f_barrierHorizontalOn.setGroup("Chemo Barrier");
  f_barrierRingOn.setGroup("Chemo Barrier");
  f_ringInnerRadius.setGroup("Chemo Barrier");
  
  f_q1Chemical.setGroup("Chemo Quadrants");
  f_q2Chemical.setGroup("Chemo Quadrants");
  f_q3Chemical.setGroup("Chemo Quadrants");
  f_q4Chemical.setGroup("Chemo Quadrants");
  f_q1Concentration.setGroup("Chemo Quadrants");
  f_q2Concentration.setGroup("Chemo Quadrants");
  f_q3Concentration.setGroup("Chemo Quadrants");
  f_q4Concentration.setGroup("Chemo Quadrants");
  
  f_staticPoints.setGroup("Chemo Static Points");
  f_staticChemicals.setGroup("Chemo Static Points");
  f_staticConcentrations.setGroup("Chemo Static Points");

  f_dynamicPoints.setGroup("Chemo Dynamic Drop Test");
  f_dynamicChemicals.setGroup("Chemo Dynamic Drop Test");
  f_dynamicConcentrations.setGroup("Chemo Dynamic Drop Test");
  f_dynamicTiming.setGroup("Chemo Dynamic Drop Test");
  
  f_dropChemicals.setGroup("Chemo Drop Test");
  f_dropVolumes.setGroup("Chemo Drop Test");
  f_dropConcentrations.setGroup("Chemo Drop Test");
  f_dropTiming.setGroup("Chemo Drop Test");

  /* .. cannot be done in constructor - the f_fileChemicalProp.getValue is empty ..
  // load chemical properties
  COUT("loadChemicalProperties: " << f_fileChemicalProp.getValue());
  //load the file that contains the physical parameters of the chemicals used
  if (!f_fileChemicalProp.getValue().empty())
  {
    m_chemicalProp = loadChemicalPorperties(f_fileChemicalProp.getFullPath().c_str());
  }  
  
  helper::set<std::string> listChemical;
  for (t_mapChemicalProp::iterator it = m_chemicalProp.begin(); it != m_chemicalProp.end(); ++it)
  {
    std::string chem = it->second.name + std::string(" ") + it->second.desc;
    listChemical.insert( chem.c_str() );
  } 

  helper::OptionsGroup chemOptions(listChemical);*/

  helper::OptionsGroup chemOptions(11,"None",
           "NaCl (attractant, water soluble)",
           "biotin (attractant, water soluble)",
           "ethanol (attractant, odorant/volatile)",
           "butanone (attractant, odorant/volatile)",
           "CuSO4 (barrier, repellent, water soluble)",
           "SDS - Sodium dodecyl sulfate (barrier, repellent, water soluble)",
           "quinine (repellent, water soluble)",
           "benzaldehyde (repellent, odorant/volatile)",
           "diacetyl (repellent, odorant/volatile)",
           "NaN3 (immobilization agent)");
  chemOptions.setSelectedItem(0);

  f_barrierChemical.setValue(chemOptions);
  f_q1Chemical.setValue(chemOptions);
  f_q2Chemical.setValue(chemOptions);
  f_q3Chemical.setValue(chemOptions);
  f_q4Chemical.setValue(chemOptions);
  
  m_isChemotaxisOn = false;
}

StimuliChemo::~StimuliChemo()
{
   
}

void StimuliChemo::init()
{
  DOUT("init StimuliChemo " << this->getName());
  Stimuli::init();

  COUT("loadChemicalProperties: " << f_fileChemicalProp.getValue());
  //load the file that contains the physical parameters of the chemicals used
  if (!f_fileChemicalProp.getValue().empty())
  {
    m_chemicalProp = loadChemicalPorperties(f_fileChemicalProp.getFullPath().c_str());
  }  

  COUT("loadSensoryPositions");
  //load the file that contains the positions of the sensory organs and their connections to the sensory-neurons
  if (!f_fileSensoryPos.getValue().empty())
  {
    m_sensoryPos = this->loadSensoryPositions(f_fileSensoryPos.getFullPath().c_str());
  }  

  m_staticPoints = f_staticPoints.getValue();
  m_staticChemicals = f_staticChemicals.getValue();
  m_staticConcentrations = f_staticConcentrations.getValue();

  m_dynamicPoints = f_dynamicPoints.getValue();
  m_dynamicChemicals = f_dynamicChemicals.getValue();
  m_dynamicConcentrations = f_dynamicConcentrations.getValue();
  m_dynamicTiming = f_dynamicTiming.getValue();
  
  m_dropChemicals = f_dropChemicals.getValue();
  m_dropVolumes = f_dropVolumes.getValue();
  m_dropConcentrations = f_dropConcentrations.getValue();
  m_dropTiming = f_dropTiming.getValue();

  m_isChemotaxisOn =   (f_barrierChemical.getValue().getSelectedId() + 
                        f_q1Chemical.getValue().getSelectedId() + 
                        f_q2Chemical.getValue().getSelectedId() + 
                        f_q3Chemical.getValue().getSelectedId() + 
                        f_q4Chemical.getValue().getSelectedId() > 0);
}

void StimuliChemo::reset() // resets the model
{
   resetStimuli();
}

StimuliChemo::t_mapChemicalProp StimuliChemo::loadChemicalPorperties(const char* filename)
{
  t_mapChemicalProp out;

  std::ifstream file(filename);
  if (!file.good()) {
    
    serr << "Chemical properties definition file \""<< filename <<"\" not found"<< sendl;
    return out;
  }

  std::string str((std::istreambuf_iterator<char>(file)),
  std::istreambuf_iterator<char>());

  Document doc;
  doc.Parse(str.c_str());
  static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  UNUSED(kTypeNames);
  
  const Value& base= doc["chemicalsCharacteristicsJSON"];
  if (!base.IsObject()){
    serr << "Error loading base chemicalsCharacteristicsJSON node from file \""<< filename <<"\"."<< sendl;
    return out;
  }
  
  const Value& cons= base["constants"];
  if (!cons.IsObject()){
    serr << "Error loading constants node from file \""<< filename <<"\"."<< sendl;
    return out;
  }
  
  for (Value::ConstMemberIterator itr = cons.MemberBegin(); itr != cons.MemberEnd(); ++itr)
  {
    //printf("Member [%s] - type is [%s] \n", itr->name.GetString(), kTypeNames[itr->value.GetType()]);
    const Value& cons_val = itr->value;

    t_chemicalProp chem;
    chem.stimulusID = cons_val["stimulusID"].GetInt();
    chem.name = cons_val["name"].GetString();
    chem.desc = cons_val["desc"].GetString();
    chem.diffuse = cons_val["diffusivity"].GetFloat();
    chem.diffuse_dt = cons_val["time_diffusion"].GetFloat();
    chem.evap_dt = cons_val["time_evaporation"].GetFloat();
    
    out[atoi(itr->name.GetString())] = chem;  
    IOUT( "ID " << atoi(itr->name.GetString()) << " stimulus " << chem.stimulusID << " name " << chem.name << " diffuse " << chem.diffuse);
  }  
  return out;
}


float StimuliChemo::getConcentrationStatic2(double r2, float c0, float lambda)// distance^2, steady concentration in the middle, diffusivity i.e. sigma
{
    // static: non-normalized gaussian used c = c0 * exp(- (r^2)/(2*lambda^2))    
    double tmp = r2/(lambda*lambda);
    return c0 * std::exp(- 0.5*tmp);
}
float StimuliChemo::getConcentrationDynamic2(double r2, float c0, float lambda, float beta, float gamma, float time)// distance^2, steady concentration in the middle, diffusivity i.e. sigma, diffusion term, evaporation term, time since start of dynamic chem point
// gamma < 1 expected
{
    // dynamic: c = c0* exp(-gamma(t1-t0)) * exp(- (r^2)/(2*lambda^2)) * exp(- (( (time)^2)/(2*beta^2))
    if (time < 0) return 0.0;    
    double tmp = r2/(lambda*lambda);
    double tmp2 = (time)/beta;
    return c0 * std::exp(-gamma*(time)) * std::exp(- 0.5*tmp) * std::exp(- 0.5*tmp2*tmp2);
}

void StimuliChemo::addChemicalSensation(std::string neuronName, int chem_type, float chem_value)
{
  if (chem_type > 0) {
    chem_type--;// modify to account for 0 indexing of chemicals in the output array, thus forgetting about the "None"
    t_chemoSensoryOutput::iterator outIt= m_chemoSensoryOutput.find(neuronName);  
    
    // at each time step the values are set to negatives (see t_chemoOutputValues.reset())
    if (outIt == m_chemoSensoryOutput.end()) {
      // if the chemical was not sensed in previous step, the new value is included in the output 
      t_chemoOutputValues outVal;
      outVal[chem_type] = chem_value;
      m_chemoSensoryOutput[neuronName] = outVal;
    }else{
    if (outIt->second[chem_type] < 0)
        outIt->second[chem_type] = chem_value; // no concentration from other scene definitions given so far, use current concentration
    else
        outIt->second[chem_type] += chem_value; // add concentrations
    } 
  }
}

void StimuliChemo::step(){
  //sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  //sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());
 
     
    DOUT("StimuliChemo compute stimuli "  << this->getName());
    for (t_mapSensoryPos::iterator it=m_sensoryPos.begin(); it != m_sensoryPos.end(); ++it)
    {
      std::string name = it->first;
      LRApos _pos = it->second;
      defaulttype::Vec3dTypes::Coord pos = WormRelPosToWorldPos(_pos.wX, _pos.wY, _pos.wZ);
 
      if (m_isChemotaxisOn){
        float chem_value = -1.0;
        int chem_type = 0;    
        double bw = f_barrierWidth.getValue()/2.0;
        
        //DOUT("StimuliChemo compute chemotaxis");  
        if (m_writeTraces) std::cout << "(CHQ)" << std::flush;
        // chemotaxis quadrants      
        if ((pos.x() > bw) && (pos.z() > bw)) // quadrant 1
        { 
          chem_value = f_q1Concentration.getValue();
          chem_type = f_q1Chemical.getValue().getSelectedId();
        } 
        else if ((pos.x() < -bw) && (pos.z() > bw)) // quadrant 2
        { 
          chem_value = f_q2Concentration.getValue();
          chem_type = f_q2Chemical.getValue().getSelectedId();
        } 
        else if ((pos.x() < -bw) && (pos.z() < -bw)) // quadrant 3
        { 
          chem_value = f_q3Concentration.getValue();
          chem_type = f_q3Chemical.getValue().getSelectedId();
        } 
        else if ((pos.x() > bw) && (pos.z() < -bw)) // quadrant 4
        { 
          chem_value = f_q4Concentration.getValue();
          chem_type = f_q4Chemical.getValue().getSelectedId();
        }
        else if ((pos.x() < bw) && (pos.x() > -bw))// barrier vertical
        {           
          chem_value = f_barrierConcentration.getValue();
          chem_type = f_barrierChemical.getValue().getSelectedId();
        }
        else if ((pos.z() < bw) && (pos.z() > -bw))// barrier horizontal
        {           
          chem_value = f_barrierConcentration.getValue();
          chem_type = f_barrierChemical.getValue().getSelectedId();
        }
		else
        {           
          chem_value = 0;
          chem_type = 0;
        }


      
        // osmotic ring
        if (f_barrierRingOn.getValue() == true){
          //DOUT("StimuliChemo  compute osmotic ring");
          if (m_writeTraces) std::cout << "(OR)" << std::flush;

          float relRad = sqrt(pos.x()*pos.x() + pos.z()*pos.z());
          relRad -= f_ringInnerRadius.getValue();
	  //std::cout << "relRad " << relRad << " bw " << bw*2 << std::endl;

          if ((relRad >= 0) && (relRad <=  bw*2.0)) {// barrier ring
            chem_value = f_barrierConcentration.getValue();
            chem_type = f_barrierChemical.getValue().getSelectedId();
          }           
        }
        
	if ((name =="ADFL") || (name=="ADFR") || (name=="PHAR") || (name=="PHAL")){
          DOUT("\n === " << name << ": "<< pos[0] << " " << pos[1] << " " << pos[2] << " -> " << " chemical " << chem_type << " value " << chem_value);
	          std::cout << "\n === " << name << ": "<< pos[0] << " " << pos[1] << " " << pos[2] << " -> " << " chemical " << chem_type << " value " << chem_value;
	}
        addChemicalSensation(name, chem_type, chem_value);
      }// end chemotaxis on
    

      if (!m_staticPoints.empty()){
        //DOUT("StimuliChemo  compute static point sources");
        if (m_writeTraces) std::cout << "(SP)" << std::flush;
      }
      for (unsigned int i=0; i< m_staticPoints.size(); i++)
      {
        int chem_type = m_staticChemicals[i];

        DOUT( "x: " << pos[0] << " y: " << pos[1] << " z: " << pos[2] );
        double distX = pos[0] - m_staticPoints[i][0];
        double distZ = pos[2] - m_staticPoints[i][1];
        double dist = distX*distX + distZ*distZ;
        
        float chem_value = getConcentrationStatic2( dist, m_staticConcentrations[i], m_chemicalProp[chem_type].diffuse);
        
        DOUT("static chemical " << chem_type << " conc " << m_staticConcentrations[i] << " diff " << m_chemicalProp[chem_type].diffuse << " value " << chem_value);
        addChemicalSensation(name, chem_type, chem_value);
      }

      if (!m_dynamicPoints.empty()){
        //DOUT("StimuliChemo  compute dynamic drop test");
        if (m_writeTraces) std::cout << "(DP)" << std::flush;
      }
      for (unsigned int i=0; i< m_dynamicPoints.size(); i++)
      {
        //this->getContext()->getTime() .. current time
        //this->getContext()->getDt() .. time step
        // TODO check time of event comes in ms while time in sofa is in seconds
	double time =  this->getContext()->getTime() - m_dynamicTiming[i]/1000;
	if (time > 0){
                 int chem_type = m_dynamicChemicals[i];

                 DOUT( "x: " << pos[0] << " y: " << pos[1] << " z: " << pos[2] );
                 double distX = pos[0] - m_dynamicPoints[i][0];
                 double distZ = pos[2] - m_dynamicPoints[i][1];
                 double dist = distX*distX + distZ*distZ;
        
                 float chem_value = getConcentrationDynamic2( dist, m_dynamicConcentrations[i], 
                       m_chemicalProp[chem_type].diffuse, m_chemicalProp[chem_type].diffuse_dt, m_chemicalProp[chem_type].evap_dt,
                       time);

	        DOUT("dynamic chemical " << chem_type << " value " << chem_value);
        	addChemicalSensation(name, chem_type, chem_value);
	}
      }
      
    }

      if (!m_dropChemicals.empty()){
        //DOUT("StimuliChemo  compute drop point test");
        if (m_writeTraces) std::cout << "(Drop)" << std::flush;
      }
      for (unsigned int i=0; i< m_dropChemicals.size(); i++)//considering dropping the chemical on the Worm's body
      {
        int chem_type = m_dropChemicals[i];
        float chem_volume = m_dropVolumes[i];
	
	double time =  this->getContext()->getTime() - m_dropTiming[i]/1000;
	if (time > 0)
        {// drop fell over the worm
	  chem_volume -= m_chemicalProp[chem_type].evap_dt * time;
	  // TODO test if in contact with the drop
          // double radius = 0.62035 * exp(chem_volume,1.0/3.0);
  
	  if (chem_volume > 0){
	        double chem_value =  m_dropConcentrations[i];
            	DOUT("chemical drop " << chem_type << " value " << chem_value);

    	        for (t_mapSensoryPos::iterator it=m_sensoryPos.begin(); it != m_sensoryPos.end(); ++it)
                {
      	 	   std::string name = it->first;
		   addChemicalSensation(name, chem_type, chem_value);
		}
	  }
       }
      }
      
    DOUT("StimuliChemo compute stimuli end " << this->getName());
}

void StimuliChemo::updatePosition(double /*dt*/)
{
}

void StimuliChemo::draw(const core::visual::VisualParams* /*vparams*/)
{
}

void StimuliChemo::exportOBJ(std::string /*name*/, std::ostream* /*out*/, std::ostream* /*mtl*/, int& /*vindex*/, int& /*nindex*/, int& /*tindex*/)
{
}

void StimuliChemo::updateVisual()
{
}

void StimuliChemo::computeBBox(const core::ExecParams*  /*params*/ )
{
}

void StimuliChemo::createStimuliString(){
  Stimuli::createStimuliString();

  DOUT("createStimuliString() - chemo " << this->getName());
  std::string chemoOut;
  for (t_chemoSensoryOutput::iterator it=m_chemoSensoryOutput.begin(); it != m_chemoSensoryOutput.end(); ++it)
  { 
    for (int i=0; i<10; i++)
      // FIXME this is currently used to lower down the network traffic, considering that the chemicals never evaporate completely,
      // and thus the concentration will never drop to 0 (if chemical is present) // TODO this does not hold for the drop case
      if (it->second[i] > 0.0) 
      {
        std::stringstream out; 
        out << "#" << m_neuronIDmap[it->first] << "," << (65535 - i) << "," << it->second[i];     
        chemoOut = chemoOut + out.str();
      }
  }
  this->m_stimuliString = chemoOut;
  DOUT("--- output stimuli str:" << this->m_stimuliString);
}

void StimuliChemo::resetStimuli(){
  
  DOUT("resetStimuli() - chemo " << this->getName());
  for (t_chemoSensoryOutput::iterator it=m_chemoSensoryOutput.begin(); it != m_chemoSensoryOutput.end(); ++it)
  {
    it->second.reset();
  }
}

} // namespace behaviormodel

} // namespace component

} // namespace sofa


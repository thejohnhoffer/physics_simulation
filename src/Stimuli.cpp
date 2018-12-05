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
#include "Stimuli.h"
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

SOFA_DECL_CLASS(Stimuli)

bool Stimuli::m_writeTracesInit = false;
bool Stimuli::initialized = false;
float Stimuli::m_plateSize = 220.0f;
std::map<std::string, int>  Stimuli::m_neuronIDmap; //name neuron , ID neuron
//std::vector<std::string> Stimuli::m_neuronList;

defaulttype::Vec3dTypes::VecCoord Stimuli::m_initialCuticle;
std::vector<unsigned int> Stimuli::m_initialCuticleIndices;
std::vector<double> Stimuli::m_minMaxXYZ(6,0.0); // storing [minX, maxX, minY, maxY, minZ, maxZ]

int StimuliClass = core::RegisterObject("Stimuli base class")
        .add< Stimuli >()
        .addLicense("LGPL")
        .addAuthor("PL")
        ;

Stimuli::Stimuli():
   f_fileNeuronIDmap(initData(&f_fileNeuronIDmap,"file_NeuronIDmap","file with Neuron to ID mapping")),
   f_plateSize( initData(&f_plateSize, (float) 66.0, "plateSize","size of the experimental Petri dish") ),
   f_writeTraces( initData(&f_writeTraces, false, "write_traces","write traces in console") ),
   f_writeTracesInit( initData(&f_writeTracesInit, false, "write_tracesInit","write initialization traces in console") )
{
   f_fileNeuronIDmap.setGroup("config files");
}

Stimuli::~Stimuli()
{
   
}

void Stimuli::init()
{
  m_writeTraces = f_writeTraces.getValue();
  DOUT("init Stimuli "  << this->getName());
  if (initialized){
    DOUT("Stimuli already initialized by previous components");
    return;
  }
  m_writeTracesInit = f_writeTracesInit.getValue();

  COUT("loadNeuronIDmap");
  //load the file that contains the mapping from the neuron name to its ID
  if (!f_fileNeuronIDmap.getValue().empty())
  {
    this->loadNeuronIDmap(f_fileNeuronIDmap.getFullPath().c_str());
  }
    
  DOUT("data model pre-computation");
  //sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  //sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());
  DOUT("cuticleOrdering");
  m_initialCuticle = cuticleOrdered();
  m_minMaxXYZ = getMinMaxCoords();
  DOUT("MinMaxCoords: " << m_minMaxXYZ[0] << " " << m_minMaxXYZ[1] << " " << m_minMaxXYZ[2] << " " << m_minMaxXYZ[3] << " "<< m_minMaxXYZ[4] << " " << m_minMaxXYZ[5]);
  DOUT("Size: " << m_minMaxXYZ[1] - m_minMaxXYZ[0] << " " << m_minMaxXYZ[3] - m_minMaxXYZ[2] << " "<< m_minMaxXYZ[5] - m_minMaxXYZ[4]);
  DOUT("Stimuli initialization finished");

  m_plateSize = f_plateSize.getValue();
  initialized = true;
}

void Stimuli::reset()
{
  m_writeTraces = f_writeTraces.getValue();
}

bool Stimuli::loadNeuronIDmap(const char* filename){
  
  std::ifstream file(filename);
  if (!file.good()) {    
    serr << "NeuronID map file \""<< filename <<"\" not found"<< sendl;
    return false;
  }

  std::string str((std::istreambuf_iterator<char>(file)),    std::istreambuf_iterator<char>());

  Document doc;
  doc.Parse(str.c_str());
  //static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  
  const Value& map= doc["neuronIDmap"];
  if (!map.IsObject()){
    serr << "Error loading neuron ID map from file \""<< filename <<"\"."<< sendl;
    return false;
  }
  
  for (Value::ConstMemberIterator itr = map.MemberBegin(); itr != map.MemberEnd(); ++itr)
  {
    IOUT("Neuron [" << itr->name.GetString() << "] - ID [" << itr->value.GetInt() <<"]");
    std::string neuronName = itr->name.GetString();
    helper::trim(neuronName);
    
    m_neuronIDmap[neuronName] = itr->value.GetInt();  
  }
  return true;
}

Stimuli::t_neuronPos Stimuli::loadNeuronPositions(const char* filename){
    t_neuronPos out;
    
  std::ifstream file(filename);
  if (!file.good()) {    
    serr << "Neuron soma position config file \""<< filename <<"\" not found"<< sendl;
    return out;
  }

  std::string str((std::istreambuf_iterator<char>(file)),  std::istreambuf_iterator<char>());

  Document doc;
  doc.Parse(str.c_str());
  //static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  
  const Value& base= doc["neuronSomaPosJSON"];
  if (!base.IsObject()){
    serr << "Error loading base neuronSomaPosJSON node from file \""<< filename <<"\"."<< sendl;
    return out;
  }
  
  for (Value::ConstMemberIterator itr = base.MemberBegin(); itr != base.MemberEnd(); ++itr)
  {
    IOUT("Neuron [" << itr->name.GetString() << "] - soma position [" << itr->value.GetFloat() <<"]");
    std::string neuronName = itr->name.GetString();
    helper::trim(neuronName);

    out[neuronName] = itr->value.GetFloat();  
  }
  return out;    
}

Stimuli::t_mapSensoryPos Stimuli::loadSensoryPositions(const char* filename){
  t_mapSensoryPos out;

  std::ifstream file(filename);
  if (!file.good()) {    
    serr << "Sensory definition file \""<< filename <<"\" not found"<< sendl;
    return out;
  }

  std::string str((std::istreambuf_iterator<char>(file)),  std::istreambuf_iterator<char>());

  Document doc;
  doc.Parse(str.c_str());
  static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  UNUSED(kTypeNames);
  
  const Value& base= doc["sensoryPosJSON"];
  if (!base.IsObject()){
    serr << "Error loading base sensoryPosJSON node from file \""<< filename <<"\"."<< sendl;
    return out;
  }
  
  const Value& cons= base["constants"];
  const Value& pos= base["positions"];
  if (!cons.IsObject() || !pos.IsObject()){
    serr << "Error loading sensory positions from file \""<< filename <<"\"."<< sendl;
    return out;
  }
  
  std::map<std::string, float> constantsMap;
  for (Value::ConstMemberIterator itr = cons.MemberBegin(); itr != cons.MemberEnd(); ++itr)
  {
    //printf("Member [%s] - type is [%s] - value is [%f] \n", itr->name.GetString(), kTypeNames[itr->value.GetType()], itr->value.GetFloat());
    constantsMap[itr->name.GetString()] = itr->value.GetFloat();  
  }

  for (Value::ConstMemberIterator itr = pos.MemberBegin(); itr != pos.MemberEnd(); ++itr)
  {
    std::string memberName = itr->name.GetString();
    helper::trim(memberName);
    //printf("Member [%s] - type is [%s] - array values: ", memberName.c_str(), kTypeNames[itr->value.GetType()]);
    LRApos pos;
    for (rapidjson::SizeType i = 0; i < itr->value.Size(); i++){
      if (itr->value[i].GetType() == 5){
        //printf(" %s,", itr->value[i].GetString());
        pos[i] = constantsMap[itr->value[i].GetString()];
      }
      else if (itr->value[i].GetType() == 6){
        //printf(" %f,", itr->value[i].GetFloat());
        pos[i] = itr->value[i].GetFloat();
      }
      else { 
        //printf(" --- "); 
      }
    }
    //printf("\n");
    WormPolarToWormRelPos(pos[0], pos[1], pos[2], pos.wX, pos.wY, pos.wZ);
    out[memberName] = pos;
    IOUT("Member " << memberName << " pos " << pos[0] << ", " << pos[1] <<", " <<pos[2]);
    IOUT( "worm rel: " << pos.wX << ", " << pos.wY << ", " << pos.wZ);
    if (m_writeTracesInit){
       defaulttype::Vec3dTypes::Coord wpos = WormRelPosToWorldPos(pos.wX, pos.wY, pos.wZ);
       IOUT( "World: " << wpos[0] << ", " << wpos[1] << ", " << wpos[2]);       
    }
  }
  
  return out;
}


///  loads the indices of the muscle rings, as specified in the muscleExtremesIndices.json config file
Stimuli::t_muscleExtremeIndices Stimuli::loadMuscleExtremes(const char* filename)
{
  t_muscleExtremeIndices out; 

  std::ifstream file(filename);
  if (!file.good()) {    
    serr << "Muscle exterme indices configuration file \""<< filename <<"\" not found"<< sendl;
    return out;
  }

  std::string str((std::istreambuf_iterator<char>(file)),  std::istreambuf_iterator<char>());

  Document doc;
  doc.Parse(str.c_str());
  static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  UNUSED(kTypeNames);
  
  const Value& base= doc["muscleExtremesIndicesJSON"];
  if (!base.IsObject()){
    serr << "Error loading base muscleExtremesIndicesJSON node from file \""<< filename <<"\"."<< sendl;
    return out;
  }
  
  out.resize(8); // 8 octants
  for (Value::ConstMemberIterator itr = base.MemberBegin(); itr != base.MemberEnd(); ++itr)
  {
    std::string octantName = itr->name.GetString();
    //helper::trim(octantName);
    int octant = std::stoi(octantName);
    //printf("Member [%s] - type is [%s] - array values: ", octantName.c_str(), kTypeNames[itr->value.GetType()]);

    std::vector<unsigned int> indices;
    for (rapidjson::SizeType i = 0; i < itr->value.Size(); i++){
      indices.push_back(itr->value[i].GetInt());
    }

    out[octant] = indices;    
  }
  
  //check file read
  if (m_writeTracesInit)
    for(unsigned int i = 0; i < out.size(); i++){
      for(unsigned int j = 0; j < out.at(i).size(); j++){
        std::cout << "octant: " << i << ", extreme point: " << j << " -> index: " << out.at(i).at(j) << std::endl;
      }
    }

  return out;
}
Stimuli::t_musclePositionConf Stimuli::loadMuscleConfiguration(const char* filename)
{
  t_musclePositionConf out;

  std::ifstream file(filename);
  if (!file.good()) {
    serr << "Muslce Positions Configurtation file \""<< filename <<"\" not found"<< sendl;
    return out;
  }

  std::string str((std::istreambuf_iterator<char>(file)),  std::istreambuf_iterator<char>());

  Document doc;
  doc.Parse(str.c_str());
  static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  UNUSED(kTypeNames);
  
  const Value& base= doc["muscleConfigurationJSON"];
  if (!base.IsObject()){
    serr << "Error loading base muscleConfigurationJSON node from file \""<< filename <<"\"."<< sendl;
    return out;
  }
  
  for (Value::ConstMemberIterator itr = base.MemberBegin(); itr != base.MemberEnd(); ++itr)
  {
    std::string muscleName = itr->name.GetString();
    helper::trim(muscleName);
    //printf("Muscle [%s] - type is [%s] - array values: ", muscleName.c_str(), kTypeNames[itr->value.GetType()]);

    t_muscleConf mConf; 
    mConf.octant = itr->value["oct"].GetInt();
    mConf.element = itr->value["mus"].GetInt();
    IOUT("Muscle " << muscleName << " octant " << mConf.octant << " element " << mConf.element);
    
    out[muscleName] = mConf;
  }
  return out;
}

Stimuli::t_muscleSensoryConf Stimuli::loadMuscleSensoryConfiguration(const char* filename)
{
  t_muscleSensoryConf out;

  std::ifstream file(filename);
  if (!file.good()) {
    serr << "MuslceSensory definition file \""<< filename <<"\" not found"<< sendl;
    return out;
  }

  std::string str((std::istreambuf_iterator<char>(file)),  std::istreambuf_iterator<char>());

  Document doc;
  doc.Parse(str.c_str());
  static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  UNUSED(kTypeNames);
  
  const Value& base= doc["sensoryMusclesJSON"];
  if (!base.IsObject()){
    serr << "Error loading base sensoryMusclesJSON node from file \""<< filename <<"\"."<< sendl;
      return out;
  }
  
  for (Value::ConstMemberIterator itr = base.MemberBegin(); itr != base.MemberEnd(); ++itr)
  {
    std::string neuronName = itr->name.GetString();
    helper::trim(neuronName);
    //printf("Neuron [%s] - type is [%s] - array values [%d]: \n", neuronName.c_str(), kTypeNames[itr->value.GetType()],itr->value.Size());

    muscleWeight mW;    
    for (rapidjson::SizeType i = 0; i < itr->value.Size(); i++){
      if (itr->value[i].GetType() == 5){
        //printf(" %s,", itr->value[i].GetString());
        mW.name = itr->value[i].GetString();
      }
      else if (itr->value[i].GetType() == 6){
        //printf(" %f,", itr->value[i].GetFloat());
        mW.weight = itr->value[i].GetFloat();
      }
      else { 
        //printf(" --- "); 
      }
    }
    //printf("\n");
    
    t_muscleSensoryConf::iterator it = out.find(neuronName);
    if (it == out.end())
      out[neuronName] = t_muscleWeightVec(1,mW);
    else 
       out[neuronName].push_back(mW);
    IOUT("neuron " << neuronName << " muscle sens size " << out[neuronName].size());
  }
  
  return out;
}


Stimuli::t_nodeNeuronDistances Stimuli::loadNodeNeuronDistances(const char* filename)
{
  /*std::string meshFilename(filename);
  if (!sofa::helper::system::DataRepository.findFile (meshFilename))
  {
    serr << "Mesh \""<< filename <<"\" not found"<< sendl;
    return false;
  }
  this->f_fileMuscleExtremes.setValue( filename );*/
  t_nodeNeuronDistances out;

  std::ifstream file(filename);
  if (!file.good()) return out;

  std::string str((std::istreambuf_iterator<char>(file)),  std::istreambuf_iterator<char>());

  Document d;
  d.Parse(str.c_str());
  static const char* kTypeNames[] = { "Null", "False", "True", "Object", "Array", "String", "Number" };
  UNUSED(kTypeNames);
  
  const Value& base= d["nodeNeuronTableJSON"];
  if (!base.IsObject()){
    serr << "Error loading base nodeNeuronTableJSON node from file \""<< filename <<"\"."<< sendl;
    return out;
  }  

  for (Value::ConstMemberIterator itr = base.MemberBegin(); itr != base.MemberEnd(); ++itr)         {             
    //printf("Member [%s] - type is [%s]\n", itr->name.GetString(), kTypeNames[itr->value.GetType()]);
    if(itr->value.IsObject()==true)
    {
      std::map<std::string, double> aux;
      IOUT( "Process Object [" << itr->name.GetString() << "]" );
      std::string nodeName = itr->name.GetString();
      int nodeID = std::stoi(nodeName);

      const Value &v = itr->value;
      for (Value::ConstMemberIterator it = v.MemberBegin(); it != v.MemberEnd(); ++it)                 
      {  
        //if(itr == base.MemberBegin()){                   
          //printf("Member [%s] is %s ", it->name.GetString(), kTypeNames[it->value.GetType()]);
        //}	
        aux[it->name.GetString()]=it->value.GetDouble();
      }

      out[nodeID]=aux;
    }
    //std::cout<<"\n"<< std::endl;
  }

  if (m_writeTracesInit)
    for(unsigned int i = 0; i < out.size(); i++){
      std::cout << "node " << i << " size of distance map " << out.at(i).size() <<std::endl;
    }

  return out;
}

std::vector<double> Stimuli::getMinMaxCoords()
{
   std::vector<double> out;
   out.push_back(m_initialCuticle[0][0]);
   out.push_back(m_initialCuticle[0][0]);
   out.push_back(m_initialCuticle[0][1]);
   out.push_back(m_initialCuticle[0][1]);
   out.push_back(m_initialCuticle[0][2]);
   out.push_back(m_initialCuticle[0][2]);

   for (unsigned int i=1; i< m_initialCuticle.size(); i++){
     out[0] = m_initialCuticle[i][0] < out[0] ? m_initialCuticle[i][0] : out[0]; 
     out[1] = m_initialCuticle[i][0] > out[1] ? m_initialCuticle[i][0] : out[1]; 
     out[2] = m_initialCuticle[i][1] < out[2] ? m_initialCuticle[i][1] : out[2]; 
     out[3] = m_initialCuticle[i][1] > out[3] ? m_initialCuticle[i][1] : out[3]; 
     out[4] = m_initialCuticle[i][2] < out[4] ? m_initialCuticle[i][2] : out[4]; 
     out[5] = m_initialCuticle[i][2] > out[5] ? m_initialCuticle[i][2] : out[5]; 
   } 
   return out;
}

defaulttype::Vec3dTypes::VecCoord Stimuli::cuticleOrdered()
{
  sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());
  defaulttype::Vec3dTypes::VecCoord aux = mechanicalObject->x0.getValue();
  defaulttype::Vec3dTypes::VecCoord aux2;

  for(unsigned int i = 0; i < 330; i++){
    bool A = true;
    for(unsigned int j = 0; j < aux2.size(); j++){
      if((aux.at(i)[0] < aux2.at(j)[0]) || (aux.at(i)[0] == aux2.at(j)[0] && aux.at(i)[1] < aux2.at(j)[1]) || (aux.at(i)[0] == aux2.at(j)[0] && aux.at(i)[1] == aux2.at(j)[1] && aux.at(i)[2] < aux2.at(j)[2])){
        aux2.insert(aux2.begin()+j,aux.at(i));
        if (!initialized)
           m_initialCuticleIndices.insert(m_initialCuticleIndices.begin()+j,i);
        A = false;
        break;
      }
    }
    if(A){
      aux2.push_back(aux.at(i));
      if (!initialized)
         m_initialCuticleIndices.push_back(i);
    }
  }

  //check ordered values
  /*for(unsigned int i = 0; i < aux2.size(); i++){
    std::cout << i << ": " << aux2.at(i) << std::endl;
  }*/

  return aux2;
}


Stimuli::t_pointVec Stimuli::computeCentralPoints()
{
  sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());
  defaulttype::Vec3dTypes::VecCoord x = mechanicalObject->x.getValue();
  
  std::vector<defaulttype::Vec3dTypes::Coord> aux;
  defaulttype::Vec3dTypes::Coord aux2;
  aux2[0] = (x[9][0] + x[4][0])/2;
  aux2[1] = (x[9][1] + x[4][1])/2;
  aux2[2] = (x[9][2] + x[4][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[308][0] + x[303][0])/2;
  aux2[1] = (x[308][1] + x[303][1])/2;
  aux2[2] = (x[308][2] + x[303][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[13][0] + x[20][0])/2;
  aux2[1] = (x[13][1] + x[20][1])/2;
  aux2[2] = (x[13][2] + x[20][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[12][0] + x[21][0])/2;
  aux2[1] = (x[12][1] + x[21][1])/2;
  aux2[2] = (x[12][2] + x[21][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[33][0] + x[40][0])/2;
  aux2[1] = (x[33][1] + x[40][1])/2;
  aux2[2] = (x[33][2] + x[40][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[32][0] + x[41][0])/2;
  aux2[1] = (x[32][1] + x[41][1])/2;
  aux2[2] = (x[32][2] + x[41][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[53][0] + x[60][0])/2;
  aux2[1] = (x[53][1] + x[60][1])/2;
  aux2[2] = (x[53][2] + x[60][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[52][0] + x[61][0])/2;
  aux2[1] = (x[52][1] + x[61][1])/2;
  aux2[2] = (x[52][2] + x[61][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[73][0] + x[80][0])/2;
  aux2[1] = (x[73][1] + x[80][1])/2;
  aux2[2] = (x[73][2] + x[80][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[72][0] + x[81][0])/2;
  aux2[1] = (x[72][1] + x[81][1])/2;
  aux2[2] = (x[72][2] + x[81][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[93][0] + x[100][0])/2;
  aux2[1] = (x[93][1] + x[100][1])/2;
  aux2[2] = (x[93][2] + x[100][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[92][0] + x[101][0])/2;
  aux2[1] = (x[92][1] + x[101][1])/2;
  aux2[2] = (x[92][2] + x[101][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[113][0] + x[120][0])/2;
  aux2[1] = (x[113][1] + x[120][1])/2;
  aux2[2] = (x[113][2] + x[120][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[112][0] + x[121][0])/2;
  aux2[1] = (x[112][1] + x[121][1])/2;
  aux2[2] = (x[112][2] + x[121][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[133][0] + x[140][0])/2;
  aux2[1] = (x[133][1] + x[140][1])/2;
  aux2[2] = (x[133][2] + x[140][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[132][0] + x[141][0])/2;
  aux2[1] = (x[132][1] + x[141][1])/2;
  aux2[2] = (x[132][2] + x[141][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[153][0] + x[160][0])/2;
  aux2[1] = (x[153][1] + x[160][1])/2;
  aux2[2] = (x[153][2] + x[160][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[152][0] + x[161][0])/2;
  aux2[1] = (x[152][1] + x[161][1])/2;
  aux2[2] = (x[152][2] + x[161][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[173][0] + x[180][0])/2;
  aux2[1] = (x[173][1] + x[180][1])/2;
  aux2[2] = (x[173][2] + x[180][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[172][0] + x[181][0])/2;
  aux2[1] = (x[172][1] + x[181][1])/2;
  aux2[2] = (x[172][2] + x[181][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[193][0] + x[200][0])/2;
  aux2[1] = (x[193][1] + x[200][1])/2;
  aux2[2] = (x[193][2] + x[200][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[192][0] + x[201][0])/2;
  aux2[1] = (x[192][1] + x[201][1])/2;
  aux2[2] = (x[192][2] + x[201][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[213][0] + x[220][0])/2;
  aux2[1] = (x[213][1] + x[220][1])/2;
  aux2[2] = (x[213][2] + x[220][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[212][0] + x[221][0])/2;
  aux2[1] = (x[212][1] + x[221][1])/2;
  aux2[2] = (x[212][2] + x[221][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[233][0] + x[240][0])/2;
  aux2[1] = (x[233][1] + x[240][1])/2;
  aux2[2] = (x[233][2] + x[240][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[232][0] + x[241][0])/2;
  aux2[1] = (x[232][1] + x[241][1])/2;
  aux2[2] = (x[232][2] + x[241][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[253][0] + x[260][0])/2;
  aux2[1] = (x[253][1] + x[260][1])/2;
  aux2[2] = (x[253][2] + x[260][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[252][0] + x[261][0])/2;
  aux2[1] = (x[252][1] + x[261][1])/2;
  aux2[2] = (x[252][2] + x[261][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[270][0] + x[280][0])/2;
  aux2[1] = (x[270][1] + x[280][1])/2;
  aux2[2] = (x[270][2] + x[280][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[273][0] + x[281][0])/2;
  aux2[1] = (x[273][1] + x[281][1])/2;
  aux2[2] = (x[273][2] + x[281][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[310][0] + x[320][0])/2;
  aux2[1] = (x[310][1] + x[320][1])/2;
  aux2[2] = (x[310][2] + x[320][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[313][0] + x[321][0])/2;
  aux2[1] = (x[313][1] + x[321][1])/2;
  aux2[2] = (x[313][2] + x[321][2])/2;
  aux.push_back(aux2);

  aux2[0] = (x[291][0] + x[296][0])/2;
  aux2[1] = (x[291][1] + x[296][1])/2;
  aux2[2] = (x[291][2] + x[296][2])/2;
  aux.push_back(aux2);

  //check 
  /*for(unsigned int i = 0; i < aux.size(); i++){
    std::cout << i << ": " << aux.at(i) << std::endl;
  }*/

  return aux;

}

// creates list of extremal muscle points (centers of the extremal muscle rings), 
// following the definition of endpoints in the muscle extreme indices (file f_fileMuscleExtremes)
Stimuli::t_muscleEnds Stimuli::computeMuscleEnds(const Stimuli::t_muscleExtremeIndices& extremeIndices)
{
  sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());
  defaulttype::Vec3dTypes::VecCoord x = mechanicalObject->x.getValue();

  //muscles ordered in rows: 0=DMR,      1=DLR,     2=VLR,     3=VMR,     4=VML,     5=VLL,     6=DLL,     7=DML 
  //muscles ordered in rows:  DMR (12), DLR (12), VLR (12), VMR (12), VML (12), VLL (11), DLL (12), DML (12)
  std::vector<std::vector<std::vector<defaulttype::Vec3dTypes::Coord> > > auxAll;

  std::vector<std::vector<defaulttype::Vec3dTypes::Coord> > auxOct;
  for(unsigned int i = 0; i < extremeIndices.size(); i++){
    auxOct.clear();
    for(unsigned int j = 1; j < extremeIndices.at(i).size() - 2; j += 2){
      std::vector<defaulttype::Vec3dTypes::Coord> auxEndpoints;// 2 endpoints of a given muscle

      defaulttype::Vec3dTypes::Coord aux3;
      aux3[0] = (x[extremeIndices.at(i).at(j)][0] + x[extremeIndices.at(i).at(j-1)][0])/2;
      aux3[1] = (x[extremeIndices.at(i).at(j)][1] + x[extremeIndices.at(i).at(j-1)][1])/2;
      aux3[2] = (x[extremeIndices.at(i).at(j)][2] + x[extremeIndices.at(i).at(j-1)][2])/2;
      auxEndpoints.push_back(aux3);

      aux3[0] = (x[extremeIndices.at(i).at(j+1)][0] + x[extremeIndices.at(i).at(j+2)][0])/2;
      aux3[1] = (x[extremeIndices.at(i).at(j+1)][1] + x[extremeIndices.at(i).at(j+2)][1])/2;
      aux3[2] = (x[extremeIndices.at(i).at(j+1)][2] + x[extremeIndices.at(i).at(j+2)][2])/2;
      auxEndpoints.push_back(aux3);

      auxOct.push_back(auxEndpoints);
    }
      auxAll.push_back(auxOct);
  }

  //check muscle extreme list
  /*for(unsigned int i = 0; i < aux.size(); i++){
    for(unsigned int j = 0; j < aux.at(i).size(); j++){
      std::cout << i << "," << j << ": " << aux.at(i).at(j) << std::endl;
    }
  }*/
  return auxAll;
}

void Stimuli::WormPolarToWormRelPos(double length, double radius, double angleDeg, double& weightX, double& weightY, double& weightZ)
{
  weightX = length;
  double angleRadians = angleDeg * M_PI/180.0;
  weightY = 0.5 + radius*sin(angleRadians)/2.0;
  weightZ = 0.5 + radius*cos(angleRadians)/2.0;
}
defaulttype::Vec3dTypes::Coord Stimuli::WormPolarToWormRelPos(double length, double radius, double angleDeg)
{
  defaulttype::Vec3dTypes::Coord out;
  WormPolarToWormRelPos(length,radius,angleDeg, out[0], out[1], out[2]);
  return out;
}

/// computes actual world coordinates from relative <0,1> coordinates in the worm frame proportion of length, height, width) .. all relative to cuticle sizes
/// !note, this implementation is specific to the worm model used
defaulttype::Vec3dTypes::Coord Stimuli::WormRelPosToWorldPos(double weightX, double weightY, double weightZ)
{
  if(weightX > 1) COUT("----- WARNING: x is not inside the worm" << weightX);
  if(weightY > 1) COUT("----- WARNING: y is not inside the worm" << weightY);
  if(weightZ > 1) COUT("----- WARNING: z is not inside the worm" << weightZ);

  // from relative X coordinates get the two nearest rings of the worm's mesh
  const int ringSize = 10; // number of nodes composing one ring of the cuticle - worm's model specific
  int nrRings = m_initialCuticle.size() / ringSize;
  int nrSegments = nrRings -1;
  int index = floor(weightX * nrSegments);
  if (index > nrSegments -1) index = nrSegments - 1; // only if length >= 1   // remember 0 indexing
  if (index < 0) index = 0; // only if length < 0
  // get ring relative X weight
  double wSegX = weightX * nrSegments - index;// iweight relative to segment = weight * nrSegments - nr segments before the index ring

  index *= ringSize;// segment to initial ring index
  //IOUT( m_initialCuticle.size() << " nrSeg " << nrSegments << " ---------------- len " << weightX << " segWX " << wSegX << " ind " << index );
  
  sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
  sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());
  defaulttype::Vec3dTypes::VecCoord currentPosition = mechanicalObject->x.getValue();

  // get representative points of the 2 rings (4 of each ring: bottom 0, top 9, left 3, right 4) - specific to worm model
  defaulttype::Vec3dTypes::Coord pr1 = currentPosition.at(m_initialCuticleIndices.at(index + 0));
  defaulttype::Vec3dTypes::Coord pr2 = currentPosition.at(m_initialCuticleIndices.at(index + 9));
  defaulttype::Vec3dTypes::Coord pr3 = currentPosition.at(m_initialCuticleIndices.at(index + 3));
  defaulttype::Vec3dTypes::Coord pr4 = currentPosition.at(m_initialCuticleIndices.at(index + 4));
  defaulttype::Vec3dTypes::Coord prr1 = currentPosition.at(m_initialCuticleIndices.at(index + 10));
  defaulttype::Vec3dTypes::Coord prr2 = currentPosition.at(m_initialCuticleIndices.at(index + 19));
  defaulttype::Vec3dTypes::Coord prr3 = currentPosition.at(m_initialCuticleIndices.at(index + 13));
  defaulttype::Vec3dTypes::Coord prr4 = currentPosition.at(m_initialCuticleIndices.at(index + 14));
  //IOUT( "pr1 : " << pr1[0] << ", " << pr1[1] << ", " << pr1[2]);
  //IOUT( "prr1: " << prr1[0] << ", " << prr1[1] << ", " << prr1[2]);
  //COUT( "w : " << wSegX << ", " << weightY << ", " << weightZ);
  
  // get corresponding Y,Z points in both rings
  defaulttype::Vec3dTypes::Coord r1 = (1-weightY)*pr1 + weightY*pr2 + (weightZ - 0.5)*(pr4-pr3) ;
  defaulttype::Vec3dTypes::Coord r2 = (1-weightY)*prr1 + weightY*prr2 + (weightZ - 0.5)*(prr4-prr3) ;
  //COUT( "r1 : " << r1[0] << ", " << r1[1] << ", " << r1[2]);
  //COUT( "r2 : " << r2[0] << ", " << r2[1] << ", " << r2[2]);
  
  // get actual position by interpolating the relative X coordinate
  defaulttype::Vec3dTypes::Coord out = (1-wSegX)*r1 + wSegX*r2;
  //COUT( "out: " << out[0] << ", " << out[1] << ", " << out[2]);
  return out; 
}

/// from worm relative polar l,r,a coordinates to the current world frame
defaulttype::Vec3dTypes::Coord Stimuli::WormPolarToWorldPos(double length, double radius, double angleDeg)
{
  defaulttype::Vec3dTypes::Coord out;
  WormPolarToWormRelPos(length,radius,angleDeg, out[0], out[1], out[2]);
  return WormRelPosToWorldPos(out[0], out[1], out[2]);
}

/// from real x,y,z coordinates in the worm frame (with worm centered according to m_minMaxXYZ) to the current world frame
// TODO not checked
defaulttype::Vec3dTypes::Coord Stimuli::WormPosToWorldPos(double x, double y, double z)
{

  defaulttype::Vec3dTypes::Coord aux;
  double weightX, weightY, weightZ;

  unsigned int i;  
  for(i = 0; i < m_initialCuticle.size()-30; i += 10){
  COUT( "i " << i);
    if(x <= m_initialCuticle.at(i+10)[0])
      break;
  }
  unsigned int index1 = i;
  unsigned int index2 = i + 19;
  unsigned int index3 = i + 5;
  unsigned int index4 = i + 14;
  COUT( "i " << i);

  //DOUT( "ActualPositionXYZ sizes = " << m_initialCuticle.size() << ", " << m_initialCuticleIndices.size());
  //DOUT( "ActualPositionXYZ indices = " << index1 << ", " << index2 << ", " << index3 << ", " << index4);

  weightX = (x - m_initialCuticle.at(index1)[0])/(m_initialCuticle.at(index2)[0] - m_initialCuticle.at(index1)[0]);
  weightY = (y - m_initialCuticle.at(index1)[1])/(m_initialCuticle.at(index2)[1] - m_initialCuticle.at(index1)[1]);
  weightZ = (z - m_initialCuticle.at(index3)[2])/(m_initialCuticle.at(index4)[2] - m_initialCuticle.at(index3)[2]);

  return WormRelPosToWorldPos(weightX, weightY, weightZ);
}




void Stimuli::createStimuliString()
{
  //include first timestep
  std::ostringstream strs;
  strs << int(this->getContext()->getTime()*1000);
  m_stimuliString = strs.str();
}


} // namespace behaviormodel

} // namespace component

} // namespace sofa


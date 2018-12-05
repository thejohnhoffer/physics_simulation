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
#ifndef SOFA_COMPONENT_BEHAVIORMODEL_EULERIANFLUID_Stimuli_H
#define SOFA_COMPONENT_BEHAVIORMODEL_EULERIANFLUID_Stimuli_H

#include <sofa/core/BehaviorModel.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/MarchingCubeUtility.h>
#include <sofa/helper/OptionsGroup.h>
#include <cctype>

// helper macros for info and debug cout
#define COUT(str) std::cout << str << std::endl;
//#define COUT(str) ;
#define DOUT(str) if (m_writeTraces) std::cout << str << std::endl;
//#define DOUT(str) ;
#define IOUT(str) if (m_writeTracesInit) std::cout << str << std::endl;
//#define IOUT(str) ;
#define UNUSED(x) (void)(x);

#define PI   3.14159265358979323846

namespace sofa
{

namespace helper {
  // trim from start (in place)
  static inline void ltrim(std::string &s) {
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  }

  // trim from end (in place)
  static inline void rtrim(std::string &s) {
      s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  }

  // trim from both ends (in place)
  static inline void trim(std::string &s) {
      ltrim(s);
      rtrim(s);
  }
}

namespace component
{

namespace behaviormodel
{
class Stimuli : public sofa::core::BehaviorModel
{
public:
  SOFA_CLASS(Stimuli, sofa::core::BehaviorModel);

  typedef std::map<std::string, float> t_neuronPos;  

  // data for encoding cuticle node to neuron distances:  nodeID -> [neuronName -> distance (nearest point of the 3D neuron model)]
  typedef std::map< std::string, double > t_neuronDistances;    
  typedef std::map<int, t_neuronDistances > t_nodeNeuronDistances;    

  /// data for encoding stimuli sensors position
  /// LRApos - Length (0..1), Radius (0..1), Angle (0..360)
  typedef struct LRApos {
    LRApos(): length(0.0),radius(0.0),angleDeg(0.0) {} 
    LRApos(float _length, float _radius, int _angleDeg):length(_length), radius(_radius), angleDeg(_angleDeg) {}
          
    float& operator[](int i) {
      if (i==0) return length;
      if (i==1) return radius;
      //if (i==2) 
      return angleDeg;// for the rest ...     
    }
    
    float length;
    float radius;
    float angleDeg;
    double wX,wY,wZ;
  } t_sensoryPos;
  typedef std::map<std::string, t_sensoryPos> t_mapSensoryPos; // ID neuron, sensory organ position 

  // extremal representative nodes (indices) of each muscle  element, grouped in octants and ordered by x (ascending)
  typedef std::vector<std::vector <unsigned int> > t_muscleExtremeIndices;
  
  // muscle sensory configuration
  // nodeName -> vector: [(muscleName, weight),...]
  typedef struct muscleWeight {
    muscleWeight(std::string n = "", float w = 0.0):name(n), weight(w) {}
    
    std::string name;
    float weight;
  } t_muscleWeight;
  typedef std::vector<t_muscleWeight> t_muscleWeightVec;
  typedef std::map<std::string, t_muscleWeightVec > t_muscleSensoryConf;
  
  // configuration of the muscle cells in respect to the 3D model of 8 octants and 12 elements each* (each apart from VLL, where only 11 are considered)
  // this configuration corresponds to the muscle ends returned by computeMuscleEnds()
  typedef struct muscleConf {
    muscleConf(int o=0, int e=0):octant(o), element(e) {}
      
    int octant;
    int element; // muscle cell in the 3D model
  } t_muscleConf;
  typedef std::map<std::string, t_muscleConf> t_musclePositionConf;

protected:
    Stimuli();
    virtual ~Stimuli();
    
public: // sofa specific mehtods

    virtual void init();
    virtual void reset();
        
    virtual void step() {}
    virtual void updatePosition(double dt) {UNUSED(dt);}
       
    virtual void draw(const core::visual::VisualParams* vparams) { UNUSED(vparams); }
    virtual void exportOBJ(std::string name, std::ostream* out, std::ostream* mtl, int& vindex, int& nindex, int& tindex) 
    { UNUSED(name); UNUSED(out); UNUSED(mtl); UNUSED(vindex); UNUSED(nindex); UNUSED(tindex); }

    virtual void updateVisual() {}
    virtual void computeBBox(const core::ExecParams*  params ) { UNUSED(params); }

public: // stimuli specific methods
  
   // helper to load neuron IDs
   ///loads neuron id map of a format "neuron_name"->(int)ID, from a JSON configuration file 
   bool loadNeuronIDmap(const char* filename);

   // helper to load for each neuron the soma position as defined on the anterior-posterior axis (normalized to 1 = length of the worm)
   t_neuronPos loadNeuronPositions(const char* filename);
   
   // helper to load for each sensory neuron the positions of its associated sensory cells
   t_mapSensoryPos loadSensoryPositions(const char* filename);

   // helpers to load muscle configuraiton files 
   t_muscleExtremeIndices loadMuscleExtremes(const char* filename);
   t_musclePositionConf loadMuscleConfiguration(const char* filename);
   t_muscleSensoryConf loadMuscleSensoryConfiguration(const char* filename);
      
   // helper to load the nearest distance from each cuticle node to each neuron
   t_nodeNeuronDistances loadNodeNeuronDistances(const char* filename); 
  
   /// initializes the ordered list of cuticle mesh points
   defaulttype::Vec3dTypes::VecCoord cuticleOrdered();
   // computes the minimum and maximum coordinates of the worm
   std::vector<double> getMinMaxCoords();
   /// the central points are hardcoded, based on the currently used 3D cuticle mesh of the worm
   typedef std::vector<defaulttype::Vec3dTypes::Coord> t_pointVec;
   t_pointVec computeCentralPoints();

   /// creates list of extremal muscle points (centers of the extremal muscle rings), 
   /// following the definition of endpoints in the muscle extreme indices (file f_fileMuscleExtremes)
   typedef std::vector<std::vector<std::vector<defaulttype::Vec3dTypes::Coord> > > t_muscleEnds;
   t_muscleEnds computeMuscleEnds(const t_muscleExtremeIndices& extremeIndices);

   /// transform polar Worm's coordinates to relative Worm's coordinates
   /// transforms Length, radius, angle coordinates to coordinates relateive to the worm's body, normalized to <0,1> (lenght, height, width)
   void WormPolarToWormRelPos(double length, double radius, double angleDeg, double& weightX, double& weightY, double& weightZ);
   defaulttype::Vec3dTypes::Coord WormPolarToWormRelPos(double length, double radius, double angleDeg);
   /// computes actual world coordinates from relative <0,1> coordinates in the worm frame (i.e. proportion of length, height, width .. all relative to cuticle sizes)
   defaulttype::Vec3dTypes::Coord WormRelPosToWorldPos(double weightX, double weightY, double weightZ);//updated computeActualPositionXYZRanges
   /// from worm relative polar l,r,a coordinates to the current world frame
   defaulttype::Vec3dTypes::Coord WormPolarToWorldPos(double length, double radius, double angleDeg); // updated computeActualPositionLRA
   /// from real x,y,z coordinates in the worm frame (with worm centered according to m_minMaxXYZ) to the current world frame
   /// TODO not ready
   defaulttype::Vec3dTypes::Coord WormPosToWorldPos(double x, double y, double z); //updated computeActualPositionXYZ

   virtual void resetStimuli() {}

   std::string getStimuliString() { return m_stimuliString; }
   void clearStimuliString() { m_stimuliString.clear(); }
   /// prints into m_stimuliString
   virtual void createStimuliString();

protected:

  sofa::core::objectmodel::DataFileName f_fileNeuronIDmap;    
  Data<float> f_plateSize;
  Data<bool> f_writeTraces;
  Data<bool> f_writeTracesInit;
  
  bool m_writeTraces;
  static bool m_writeTracesInit;
  static bool initialized;
  static float m_plateSize;
  static std::map<std::string, int>  m_neuronIDmap; //name neuron , ID neuron   
  
  // the following is used when translating positions from local worm coordinates to world coordinates
  static defaulttype::Vec3dTypes::VecCoord m_initialCuticle;
  static std::vector<unsigned int> m_initialCuticleIndices;

  /// storing [minX, maxX, minY, maxY, minZ, maxZ]
  /// Current model MinMaxCoords: -5.6221 6.0668 -0.092 1.2078 1.0725 2.349
  /// Current Model Size: 11.6889 1.2998 1.2765
  static std::vector<double> m_minMaxXYZ; 

  
  // formatted output string of the generated stimuli - for the base class only time is generated
  std::string m_stimuliString; 
};


} // namespace behaviormodel

} // namespace component

} // namespace sofa

#endif

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
#ifndef SOFA_COMPONENT_BEHAVIORMODEL_StimuliChemo_H
#define SOFA_COMPONENT_BEHAVIORMODEL_StimuliChemo_H

#include <sofa/core/BehaviorModel.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/MarchingCubeUtility.h>
#include <sofa/helper/OptionsGroup.h>
#include "Stimuli.h"
#include <sofa/defaulttype/Vec.h>
//#include <sofa/helper/vector.h>
//#include <sofa/defaulttype/VecTypes.h>

namespace sofa
{

namespace component
{

namespace behaviormodel
{

class StimuliChemo : public sofa::component::behaviormodel::Stimuli
{
public:
    SOFA_CLASS(StimuliChemo, sofa::component::behaviormodel::Stimuli);

protected:
public:
protected:
    StimuliChemo();
    virtual ~StimuliChemo();
    
public: // sofa specific mehtods

    virtual void init();
    virtual void reset();
        
    virtual void step();
    virtual void updatePosition(double dt);
       
    virtual void draw(const core::visual::VisualParams* vparams);
    virtual void exportOBJ(std::string name, std::ostream* out, std::ostream* mtl, int& vindex, int& nindex, int& tindex);

    virtual void updateVisual();
    virtual void computeBBox(const core::ExecParams*  params );

public: // stimuli specific methods
  
  virtual void resetStimuli();

  virtual void createStimuliString(); 

protected:
  // called from step to set the sensed value into the output matrix
      void addChemicalSensation(std::string neuronName, int chem_type, float chem_value);

  struct t_chemicalProp{
    t_chemicalProp(int _stimID = 0, std::string _name="None", std::string _desc="", float _diffuse=0.0, float _diffuse_dt = 0.0, float _evap_dt = 0.0):
                   stimulusID(_stimID), name(_name), desc(_desc), diffuse(_diffuse), diffuse_dt(_diffuse_dt), evap_dt(_evap_dt) { }
    int stimulusID; 
    std::string name; 
    std::string desc; 
    float diffuse; 
    float diffuse_dt; 
    float evap_dt;
  };
  typedef std::map<int, t_chemicalProp> t_mapChemicalProp;
  t_mapChemicalProp m_chemicalProp; 

  t_mapChemicalProp loadChemicalPorperties(const char* filename);

  // distance^2, steady concentration at the middle, diffusivity i.e. sigma
  float getConcentrationStatic2(double r2, float c0, float lambda);
  
  // distance^2, steady concentration at the middle, diffusivity i.e. sigma, diffusion term, evaporation term, time since start of dynamic chem point event
  // gamma < 1 expected
  float getConcentrationDynamic2(double r2, float c0, float lambda, float beta, float gamma, float time);

protected:
  typedef Stimuli::LRApos LRApos;
  typedef Stimuli::t_mapSensoryPos t_mapSensoryPos;
  
  sofa::core::objectmodel::DataFileName f_fileSensoryPos;
  sofa::core::objectmodel::DataFileName f_fileChemicalProp;
   
  Data<float> f_barrierWidth;// relative to worm size ?
  Data<sofa::helper::OptionsGroup> f_barrierChemical;// chemical ID
  Data<float> f_barrierConcentration;// molar concentration in molar = M = mole/litre
  Data<bool> f_barrierVerticalOn;
  Data<bool> f_barrierHorizontalOn;
  Data<bool> f_barrierRingOn;
  Data<float> f_ringInnerRadius;// relative to worm size ?
  
  Data<sofa::helper::OptionsGroup> f_q1Chemical;// chemical ID
  Data<sofa::helper::OptionsGroup> f_q2Chemical;
  Data<sofa::helper::OptionsGroup> f_q3Chemical;
  Data<sofa::helper::OptionsGroup> f_q4Chemical;
  Data<float> f_q1Concentration;// molar concentration in molar = M = mole/litre
  Data<float> f_q2Concentration;
  Data<float> f_q3Concentration;
  Data<float> f_q4Concentration;

  bool m_isChemotaxisOn;
  
  // static chemo test considering constant circular (gaussian) gradient of the chemical
  Data<helper::vector<sofa::defaulttype::Vec2d> > f_staticPoints;// 2D point in the plate coordinates
  Data<sofa::helper::vector<int> > f_staticChemicals;// chemical ID
  Data<sofa::helper::vector<double> > f_staticConcentrations;// molar concentration in molar = M = mole/litre
  
  // dynamic chemo test considering circular (gaussian) gradient, diffusion and evaporation of the chemical
  Data<helper::vector<sofa::defaulttype::Vec2d> > f_dynamicPoints;// 2D point in the plate coordinates
  Data<sofa::helper::vector<int> > f_dynamicChemicals;// chemical ID
  Data<sofa::helper::vector<double> > f_dynamicConcentrations;// molar concentration in molar = M = mole/litre
  Data<sofa::helper::vector<double> > f_dynamicTiming;// in ms
 
  // drop test applying a drop of a chemical on the worm's body 
  Data<sofa::helper::vector<int> > f_dropChemicals;// chemical ID
  Data<sofa::helper::vector<double> > f_dropVolumes;// volume in mL
  Data<sofa::helper::vector<double> > f_dropConcentrations;// molar concentration in molar = M = mole/litre
  Data<sofa::helper::vector<double> > f_dropTiming;// in ms

  helper::vector<sofa::defaulttype::Vec2d> m_staticPoints;
  sofa::helper::vector<int> m_staticChemicals;
  sofa::helper::vector<double> m_staticConcentrations;

  helper::vector<sofa::defaulttype::Vec2d> m_dynamicPoints;
  sofa::helper::vector<int> m_dynamicChemicals;
  sofa::helper::vector<double> m_dynamicConcentrations;
  sofa::helper::vector<double> m_dynamicTiming;

  sofa::helper::vector<int> m_dropChemicals;
  sofa::helper::vector<double> m_dropVolumes;
  sofa::helper::vector<double> m_dropConcentrations;
  sofa::helper::vector<double> m_dropTiming;
  
  t_mapSensoryPos m_sensoryPos; // ID neuron, sensory organ position  

  
  struct t_chemoOutputValues{
    t_chemoOutputValues(int _size = 10):size(_size) { data.resize(_size,-1.0f); }
    void reset() { data.assign(size, -1.0f); }
    float& operator[](int i) { if ((i>size) || (i<0)) return dummy; return data[i]; }
    std::vector<float> data;
    static float dummy;
    int size;
  };
  typedef std::map<std::string, t_chemoOutputValues> t_chemoSensoryOutput;
  t_chemoSensoryOutput m_chemoSensoryOutput;  
  
};


} // namespace behaviormodel

} // namespace component

} // namespace sofa

#endif

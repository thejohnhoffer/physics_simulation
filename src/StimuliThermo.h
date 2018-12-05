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
#ifndef SOFA_COMPONENT_BEHAVIORMODEL_StimuliThermo_H
#define SOFA_COMPONENT_BEHAVIORMODEL_StimuliThermo_H

#include <sofa/core/BehaviorModel.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/MarchingCubeUtility.h>
#include <sofa/helper/OptionsGroup.h>
#include "Stimuli.h"

namespace sofa
{

namespace component
{

namespace behaviormodel
{

class StimuliThermo : public sofa::component::behaviormodel::Stimuli
{
public:
    SOFA_CLASS(StimuliThermo, sofa::component::behaviormodel::Stimuli);

protected:
public:
protected:
    StimuliThermo();
    virtual ~StimuliThermo();
    
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
  // distance^2, steady temperature in the middle, diffusivity i.e. sigma
  float getTemepratureStatic2(double r2, float t0, float sigma);
   
  virtual void resetStimuli();

  virtual void createStimuliString();

protected:
  typedef Stimuli::t_pointVec t_pointVec;
  typedef Stimuli::LRApos LRApos;
  //typedef Stimuli::t_sensoryPos t_sensoryPos;
  typedef Stimuli::t_neuronPos t_neuronPos;

  sofa::core::objectmodel::DataFileName f_fileNeuronSomaPos;
  
  Data<float> f_envTemperature;
  Data<float> f_heatDiffusion;
  Data<bool>  f_gradientOn;
  Data<float>  f_LRgradDistance;
  Data<float>  f_gradLeft;
  Data<float>  f_gradRight;  
  Data<helper::vector<sofa::defaulttype::Vec2d> > f_heatPoints;
  Data<sofa::helper::vector<double> > f_heatPointTemperature;
  Data<sofa::helper::vector<double> > f_heatPointTiming;// in ms
  Data<sofa::helper::vector<double> > f_heatPointDuration;// in ms
  
  float m_envTemperature;
  float m_heatDiffusion;
  float m_LRgradDistance;
  float m_gradLeft;
  float m_gradRight;
  bool m_heatGradientOn;
  helper::vector<sofa::defaulttype::Vec2d> m_heatPoints;
  sofa::helper::vector<double> m_heatPointTemperature;  
  sofa::helper::vector<double> m_heatPointTiming;  
  sofa::helper::vector<double> m_heatPointDuration;  
  
  typedef std::map< std::string, float > t_temepratureSensoryOut;
  t_temepratureSensoryOut m_neuronTemepratureSensoryOut;
  t_neuronPos m_neuronSomaPos;

};


} // namespace behaviormodel

} // namespace component

} // namespace sofa

#endif

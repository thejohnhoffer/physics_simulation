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
#ifndef SOFA_COMPONENT_BEHAVIORMODEL_[Stimulus]_H
#define SOFA_COMPONENT_BEHAVIORMODEL_[Stimulus]_H

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

class [Stimulus] : public sofa::component::behaviormodel::Stimuli
{
public:
    SOFA_CLASS([Stimulus], sofa::component::behaviormodel::Stimuli);

protected:
public:
protected:
    [Stimulus]();
    virtual ~[Stimulus]();
    
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
};


} // namespace behaviormodel

} // namespace component

} // namespace sofa

#endif

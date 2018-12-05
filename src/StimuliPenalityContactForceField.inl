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
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_INTERACTIONFORCEFIELD_StimuliPenalityContactForceField_INL
#define SOFA_COMPONENT_INTERACTIONFORCEFIELD_StimuliPenalityContactForceField_INL

#include "StimuliPenalityContactForceField.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/system/config.h>
#include <cassert>
#include <sofa/helper/gl/template.h>
#include <iostream>

#include <sofa/simulation/common/Simulation.h>
#include <SofaBaseMechanics/MechanicalObject.h>


namespace sofa
{

namespace component
{

namespace interactionforcefield
{

template<class DataTypes>
StimuliPenalityContactForceField<DataTypes>::StimuliPenalityContactForceField(MechanicalState* object1, MechanicalState* object2)	
        : Inherit(object1, object2), contacts(initData(&contacts,"contacts", "Contacts"))
    {
		m_stimuli = NULL;
    }


template<class DataTypes>
void StimuliPenalityContactForceField<DataTypes>::clear(int reserve)
{
    prevContacts.swap(*contacts.beginEdit()); // save old contacts in prevContacts
    contacts.beginEdit()->clear();
    if (reserve)
        contacts.beginEdit()->reserve(reserve);
    contacts.endEdit();
}


template<class DataTypes>
void StimuliPenalityContactForceField<DataTypes>::addContact(int m1, int m2, int index1, int index2, const Deriv& norm, Real dist, Real ks, Real /*mu_s*/, Real /*mu_v*/, int oldIndex)
{
    int i = contacts.getValue().size();
    contacts.beginEdit()->resize(i+1);
    Contact& c = (*contacts.beginEdit())[i];
    c.m1 = m1;
    c.m2 = m2;
    c.index1 = index1;
    c.index2 = index2;
    c.norm = norm;
    c.dist = dist;
    c.ks = ks;
//	c.mu_s = mu_s;
//	c.mu_v = mu_v;
    c.pen = 0;
    if (oldIndex > 0 && oldIndex <= (int)prevContacts.size())
    {
        c.age = prevContacts[oldIndex-1].age+1;
    }
    else
    {
        c.age = 0;
    }
    contacts.endEdit();
}

template<class DataTypes>
void StimuliPenalityContactForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& data_f1, DataVecDeriv& data_f2, const DataVecCoord& data_x1, const DataVecCoord& data_x2, const DataVecDeriv& /*data_v1*/, const DataVecDeriv& /*data_v2*/ )
{
	if(m_stimuli == NULL){
		sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
		sofa::simulation::Node* parent2 = dynamic_cast<sofa::simulation::Node*> (parent->getParents().at(0));
		m_stimuli = parent2->getTreeObject<behaviormodel::StimuliMechano>();
		if (!m_stimuli) COUT("binding to StimuliMechano falied");
	}
    VecDeriv&       f1 = *data_f1.beginEdit();
    const VecCoord& x1 =  data_x1.getValue();
    //const VecDeriv& v1 =  data_v1.getValue();
    VecDeriv&       f2 = *data_f2.beginEdit();
    const VecCoord& x2 =  data_x2.getValue();
    //const VecDeriv& v2 =  data_v2.getValue();

    f1.resize(x1.size());
    f2.resize(x2.size());
    for (unsigned int i=0; i<contacts.getValue().size(); i++)
    {
        Contact& c = (*contacts.beginEdit())[i];
        Coord u = x2[c.m2]-x1[c.m1];
		
        c.pen = c.dist - u*c.norm;
        if (c.pen > 0)
        {

			sofa::simulation::Node* parent = dynamic_cast<sofa::simulation::Node*> (this->getContext());
			sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>* mechanicalObject = dynamic_cast<sofa::component::container::MechanicalObject<defaulttype::Vec3dTypes>*> (parent->getMechanicalState());
			defaulttype::Vec3dTypes::VecCoord aux = mechanicalObject->x.getValue();
			/*for(unsigned int i = 0; i < x2.size(); i++){
				std::cout << i << ": " << x2.at(i) << "a" << std::endl;
			}
			for(unsigned int i = 0; i < x1.size(); i++){
				std::cout << i << ": " << x1.at(i) << "b" << std::endl;
			}
			for(unsigned int i = 0; i < 330; i++){
				std::cout << i << ": " << aux.at(i) << "c" << std::endl;
			}*/
            Real fN = c.ks * c.pen;
            Deriv force = -c.norm*fN;

            f1[c.m1]+=force;
            this->mask1->insertEntry(c.m1);
            f2[c.m2]-=force;
            this->mask2->insertEntry(c.m2);
		if (m_stimuli){
			// printf("sending forces to StimuliMechano\n");
			if(c.index1 < c.index2 && c.index1 < 330){
				m_stimuli->addMechanoStimuli(c.index1,1000000000000*sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]));
			}
			else if(c.index2 < 330){
				m_stimuli->addMechanoStimuli(c.index2,1000000000000*sqrt(force[0]*force[0]+force[1]*force[1]+force[2]*force[2]));
			}
		}
        }
    }
    contacts.endEdit();

    data_f1.endEdit();
    data_f2.endEdit();

}

template<class DataTypes>
void StimuliPenalityContactForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& data_df1, DataVecDeriv& data_df2, const DataVecDeriv& data_dx1, const DataVecDeriv& data_dx2)
{
	//std::cout << "addDForce" << std::endl;
    VecDeriv&        df1 = *data_df1.beginEdit();
    VecDeriv&        df2 = *data_df2.beginEdit();
    const VecDeriv&  dx1 =  data_dx1.getValue();
    const VecDeriv&  dx2 =  data_dx2.getValue();
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
	//std::cout << "kfactor: " << kFactor << std::endl;
    df1.resize(dx1.size());
    df2.resize(dx2.size());
	//std::cout << "sizes: " << dx1.size() << " " << dx2.size() << std::endl;
	//std::cout << "contact size: " << contacts.getValue().size() << std::endl;
    for (unsigned int i=0; i<contacts.getValue().size(); i++)
    {
        const Contact& c = contacts.getValue()[i];
		//std::cout << "aber: " << c.pen << " " << c.norm << " " << c.ks << "\n" << c.m1 << " " << dx1[c.m1] << " " << c.m2 << " " << dx2[c.m2] << std::endl;
        if (c.pen > 0) // + dpen > 0)
        {
            Coord du = dx2[c.m2]-dx1[c.m1];
            Real dpen = - du*c.norm;
            //if (c.pen < 0) dpen += c.pen; // start penality at distance 0
            Real dfN = c.ks * dpen * (Real)kFactor;
            Deriv dforce = -c.norm*dfN;
            df1[c.m1]+=dforce;
            df2[c.m2]-=dforce;
        }
    }

    data_df1.endEdit();
    data_df2.endEdit();

}

template <class DataTypes>
SReal StimuliPenalityContactForceField<DataTypes>::getPotentialEnergy(const sofa::core::MechanicalParams*, const DataVecCoord&, const DataVecCoord& ) const
{
    serr<<"StimuliPenalityContactForceField::getPotentialEnergy-not-implemented !!!"<<sendl;
    return 0;
}

template<class DataTypes>
void StimuliPenalityContactForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!((this->mstate1 == this->mstate2)?vparams->displayFlags().getShowForceFields():vparams->displayFlags().getShowInteractionForceFields())) return;
    const VecCoord& p1 = this->mstate1->read(core::ConstVecCoordId::position())->getValue();
    const VecCoord& p2 = this->mstate2->read(core::ConstVecCoordId::position())->getValue();
    //glDisable(GL_LIGHTING); // do not use gl under draw component, it cause crash when using other render !

    std::vector< defaulttype::Vector3 > points[4];

    for (unsigned int i=0; i<contacts.getValue().size(); i++)
    {
        const Contact& c = contacts.getValue()[i];
        Real d = c.dist - (p2[c.m2]-p1[c.m1])*c.norm;
        if (c.age > 10) //c.spen > c.mu_s * c.ks * 0.99)
            if (d > 0)
            {
                points[0].push_back(p1[c.m1]);
                points[0].push_back(p2[c.m2]);
            }
            else
            {
                points[1].push_back(p1[c.m1]);
                points[1].push_back(p2[c.m2]);
            }
        else if (d > 0)
        {
            points[2].push_back(p1[c.m1]);
            points[2].push_back(p2[c.m2]);
        }
        else
        {
            points[3].push_back(p1[c.m1]);
            points[3].push_back(p2[c.m2]);
        }
    }
    vparams->drawTool()->drawLines(points[0], 1, defaulttype::Vec<4,float>(1,0,1,1));
    vparams->drawTool()->drawLines(points[1], 1, defaulttype::Vec<4,float>(0,1,1,1));
    vparams->drawTool()->drawLines(points[2], 1, defaulttype::Vec<4,float>(1,0,0,1));
    vparams->drawTool()->drawLines(points[3], 1, defaulttype::Vec<4,float>(0,1,0,1));


    std::vector< defaulttype::Vector3 > pointsN;
    if (vparams->displayFlags().getShowNormals())
    {
        for (unsigned int i=0; i<contacts.getValue().size(); i++)
        {
            const Contact& c = contacts.getValue()[i];
            Coord p = p1[c.m1] - c.norm;
            pointsN.push_back(p1[c.m1]);
            pointsN.push_back(p);


            p = p2[c.m2] + c.norm;
            pointsN.push_back(p2[c.m2]);
            pointsN.push_back(p);
        }
        vparams->drawTool()->drawLines(pointsN, 1, defaulttype::Vec<4,float>(1,1,0,1));
    }
}


template<class DataTypes>
void StimuliPenalityContactForceField<DataTypes>::grabPoint(
    const core::behavior::MechanicalState<defaulttype::Vec3Types> *tool,
    const helper::vector< unsigned int > &index,
    helper::vector< std::pair< core::objectmodel::BaseObject*, defaulttype::Vec3f> > &result,
    helper::vector< unsigned int > &triangle,
    helper::vector< unsigned int > &index_point)
{
    if (static_cast< core::objectmodel::BaseObject *>(this->mstate1) == static_cast< const core::objectmodel::BaseObject *>(tool))
    {
        for (unsigned int i=0; i<contacts.getValue().size(); i++)
        {
            for (unsigned int j=0; j<index.size(); j++)
            {
                if (contacts.getValue()[i].m1  == (int)index[j])
                {
                    result.push_back(std::make_pair(static_cast< core::objectmodel::BaseObject *>(this),
                            (this->mstate2->read(core::ConstVecCoordId::position())->getValue())[contacts.getValue()[i].m2])
                                    );
                    triangle.push_back(contacts.getValue()[i].index2);
                    index_point.push_back(index[j]);
                }
            }
        }
    }
    else if (static_cast< core::objectmodel::BaseObject *>(this->mstate2) == static_cast< const core::objectmodel::BaseObject *>(tool))
    {

        for (unsigned int i=0; i<contacts.getValue().size(); i++)
        {
            for (unsigned int j=0; j<index.size(); j++)
            {
                if (contacts.getValue()[i].m2  == (int)index[j])
                {
                    result.push_back(std::make_pair(static_cast< core::objectmodel::BaseObject *>(this),
                            (this->mstate1->read(core::ConstVecCoordId::position())->getValue())[contacts.getValue()[i].m1])
                                    );

                    triangle.push_back(contacts.getValue()[i].index1);
                    index_point.push_back(index[j]);
                }
            }
        }
    }


}


} // namespace interactionforcefield

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_INTERACTIONFORCEFIELD_StimuliPenalityContactForceField_INL */




#include "Friction.h"

#include <SofaMiscFem/BoyceAndArruda.h>
#include <SofaMiscFem/NeoHookean.h>
#include <SofaMiscFem/MooneyRivlin.h>
#include <SofaMiscFem/VerondaWestman.h>
#include <SofaMiscFem/STVenantKirchhoff.h>
#include <SofaMiscFem/HyperelasticMaterial.h>
#include <SofaMiscFem/Costa.h>
#include <SofaMiscFem/Ogden.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <sofa/helper/gl/template.h>
#ifndef SOFA_NO_OPENGL
#if defined (__APPLE__)
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#endif
#include <sofa/core/behavior/ForceField.inl>
#include <SofaBaseTopology/TopologyData.inl>
#include <algorithm>
#include <iterator>
#include <sofa/helper/AdvancedTimer.h>
#include <string>
#include "CGLinearSolverTraces.h"
#include <Eigen/src/LU/FullPivLU.h>


namespace sofa
{
	namespace component
	{
		namespace forcefield
		{
			using namespace sofa::defaulttype;
			using namespace	sofa::component::topology;
			using namespace core::topology;

#define PERIOD 1.6f

			template <class DataTypes> 
			Friction<DataTypes>::Friction()
				: m_mu(initData(&m_mu, "mu", "mu"))
				, m_floorK ( initData(&m_floorK, (Real)0.0, "floor_K", "floor_K") )
				, m_springK( initData(&m_springK, "spring_K", "mass density that allows to compute the  particles masses from a mesh topology and geometry.\nOnly used if > 0") )
				, m_floorHeight ( initData(&m_floorHeight, (Real)0.0, "floorHeight", "floorHeight") )
				, f_writeTraces( initData(&f_writeTraces, true, "write_traces","write traces in console") )
				, f_writeTracesD( initData(&f_writeTracesD, true, "write_traces_diff","write diff traces in console") )
				, f_fileFriction(initData(&f_fileFriction,"friction_file","file with relations for friction"))
			{
				m_firstStep = true;
			}

			template <class DataTypes> 
			Friction<DataTypes>::~Friction()
			{
			}



			template<class DataTypes>
			void Friction<DataTypes>::draw(const core::visual::VisualParams* vparams)
			{
			}

			template <class DataTypes> 
			void Friction<DataTypes>::init()
			{
				std::cout << "init Friction" << std::endl;
				this->Inherited::init();
				if (!f_fileFriction.getValue().empty())
				{
					this->load(f_fileFriction.getFullPath().c_str());
				}
				frogSum.push_back(0.0);
				frogSum.push_back(0.0);
				frogSum.push_back(0.0);
			} // init()

			template <class DataTypes> 
			void Friction<DataTypes>::reset()
			{
				m_firstStep = true;
				m_nodes.clear();
				frogSum.at(0) = 0.0;
				frogSum.at(1) = 0.0;
				frogSum.at(2) = 0.0;
			} // init()

			template <class DataTypes> 
			bool Friction<DataTypes>::load(const char* filename)
			{
				/*std::string meshFilename(filename);
				if (!sofa::helper::system::DataRepository.findFile (meshFilename))
				{
					serr << "Mesh \""<< filename <<"\" not found"<< sendl;
					return false;
				}
				this->f_fileFriction.setValue( filename );*/

				std::ifstream file(filename);
				if (!file.good()) return false;

				std::string line;
				

				while(std::getline(file, line)){
					std::vector<std::string> values = split(line,' ');
					if(values.size() == 2){
						m_frictionRelations[std::stoi(values.at(0))] = std::stoi(values.at(1));
					}
					else if(values.size() == 1){
						m_frictionRelations[std::stoi(values.at(0))] = -1;
					}					
				}

				//check
				for (std::map<int,int>::iterator it=m_frictionRelations.begin(); it!=m_frictionRelations.end(); ++it){
					std::cout << it->first << " " << it->second << "\n";
				}

				return true;
			}

template <class DataTypes>
std::vector<std::string> &Friction<DataTypes>::split(std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

template <class DataTypes>
std::vector<std::string> Friction<DataTypes>::split(std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

template <class DataTypes> 
void Friction<DataTypes>::addForce(const core::MechanicalParams*  mparams  /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv&  d_v )
{
	const VecCoord& mu = m_mu.getValue();
	const VecCoord& springK = m_springK.getValue();

	if(mu.size() < 1 || springK.size() < 1) return;
	bool writeTraces = f_writeTraces.getValue();
	//initialize colliding vector
	const VecCoord& x = d_x.getValue();
	//std::cout << "friction addForce " << x[228] << std::endl;
	VecDeriv& f = *d_f.beginEdit();

	std::vector<double> frog;
	frog.push_back(0);
	frog.push_back(0);
	frog.push_back(0);

	const VecCoord& v = d_v.getValue();

	std::vector<double> frictionForceDirection;
	std::vector<double> frictionForcePerpendicular;

	std::vector<double> frictionForceSum;
	frictionForceSum.push_back(0);
	frictionForceSum.push_back(0);
	frictionForceSum.push_back(0);

	if(m_firstStep){
		for(unsigned int i = 0; i < x.size(); i++){
			nodeProperties aux;
			aux.m_anchor.push_back(0);
			aux.m_anchor.push_back(0);
			aux.m_anchor.push_back(0);

			aux.m_anchorDifference.push_back(0);
			aux.m_anchorDifference.push_back(0);
			aux.m_anchorDifference.push_back(0);

			aux.m_anchorDifferenceUnit.push_back(0);
			aux.m_anchorDifferenceUnit.push_back(0);
			aux.m_anchorDifferenceUnit.push_back(0);

			aux.m_directionVector.push_back(0);
			aux.m_directionVector.push_back(0);
			aux.m_directionVector.push_back(0);

			aux.m_penalty = 0.0;
			aux.m_collisionState = NOTCOLLIDE_NOTCOLLIDE;
			m_nodes.push_back(aux);
		}
		m_firstStep = false;
	}

	for(unsigned int i = 0; i < x.size(); i++){
		if(m_frictionRelations[i] && x[i][1] <= m_floorHeight.getValue()){

			f[i][2] += 10;
			std::vector<double> auxForce;
			auxForce.push_back(0);
			auxForce.push_back(0);
			auxForce.push_back(0);

			if(i == 0 && writeTraces) std::cout << "anchor: " << m_nodes[i].m_anchor[0] << " " << m_nodes[i].m_anchor[1] << " " << m_nodes[i].m_anchor[2] << std::endl;
			if(m_nodes[i].m_collisionState > 4){
				if(i == 0 && writeTraces) std::cout << "was not Colliding" << std::endl;

				m_nodes[i].m_anchor[0] = x[i][0];
				m_nodes[i].m_anchor[1] = 0;
				m_nodes[i].m_anchor[2] = x[i][2];
				m_nodes[i].m_collisionState = NOTCOLLIDE_COLLIDE;
			}
			else{

				//penalty force
				if(i == 0 && writeTraces) std::cout << "was Colliding: " << f[i] << std::endl;
				m_nodes[i].m_penalty = -m_floorK.getValue()*(x[i][1] - m_floorHeight.getValue());
							
				auxForce[1] += m_nodes[i].m_penalty;
				if(i == 0 && writeTraces) std::cout << "penalty: " << m_nodes[i].m_penalty << std::endl;

				if(m_frictionRelations.at(i) > 0){

					//compute the direction of the body at this point
					m_nodes[i].m_directionVector[0] = x[m_frictionRelations.at(i)][0] - x[i][0];
					m_nodes[i].m_directionVector[2] = x[m_frictionRelations.at(i)][2] - x[i][2];
					m_nodes[i].m_directionNorm = sqrt(m_nodes[i].m_directionVector[0]*m_nodes[i].m_directionVector[0]+m_nodes[i].m_directionVector[2]*m_nodes[i].m_directionVector[2]);
					m_nodes[i].m_directionVector[0] = m_nodes[i].m_directionVector[0]/m_nodes[i].m_directionNorm;
					m_nodes[i].m_directionVector[2] = m_nodes[i].m_directionVector[2]/m_nodes[i].m_directionNorm;

					if(i == 0 && writeTraces) std::cout << "position: " << i << ": " << x[i][0] << " " << x[i][1] << " " << x[i][2] << std::endl;
					if(i == 0 && writeTraces) std::cout << "direction node: " << m_frictionRelations.at(i) << ": " << x[m_frictionRelations.at(i)][0] << " " << x[m_frictionRelations.at(i)][2] << std::endl;
					if(i == 0 && writeTraces) std::cout << "direction: " << m_nodes[i].m_directionVector.at(0) << " " << m_nodes[i].m_directionVector.at(2) << std::endl;
					if(i == 0 && writeTraces) std::cout << "direction module: " << m_nodes[i].m_directionNorm << std::endl;

					frictionForcePerpendicular.clear();
					frictionForceDirection.clear();

					//first difference between current point and anchor
					m_nodes[i].m_anchorDifference[0] = x[i][0] - m_nodes[i].m_anchor[0]; 
					m_nodes[i].m_anchorDifference[2] = x[i][2] - m_nodes[i].m_anchor[2];
					m_nodes[i].m_differenceNorm2 = m_nodes[i].m_anchorDifference[0]*m_nodes[i].m_anchorDifference[0] + m_nodes[i].m_anchorDifference[2]*m_nodes[i].m_anchorDifference[2];
					m_nodes[i].m_differenceNorm = sqrt(m_nodes[i].m_differenceNorm2);
					m_nodes[i].m_anchorDifferenceUnit[0] = m_nodes[i].m_anchorDifference[0]/m_nodes[i].m_differenceNorm; 
					m_nodes[i].m_anchorDifferenceUnit[2] = m_nodes[i].m_anchorDifference[2]/m_nodes[i].m_differenceNorm; 

					if(i == 0 && writeTraces) std::cout << "anchor: " << m_nodes[i].m_anchor[0] << " " << m_nodes[i].m_anchor[2] << std::endl;
					if(i == 0 && writeTraces) std::cout << "difference: " << m_nodes[i].m_anchorDifference[0] << " " << m_nodes[i].m_anchorDifference[2] << "; modules: " << m_nodes[i].m_differenceNorm << " " << m_nodes[i].m_differenceNorm2 << std::endl;
					if(i == 0 && writeTraces) std::cout << "differenceUnit: " << m_nodes[i].m_anchorDifferenceUnit[0] << " " << m_nodes[i].m_anchorDifferenceUnit[2] << std::endl;

					double scalarDirection = m_nodes[i].m_anchorDifference[0]*m_nodes[i].m_directionVector[0] + m_nodes[i].m_anchorDifference[2]*m_nodes[i].m_directionVector[2];
					double scalarPerpendicular = m_nodes[i].m_anchorDifference[0]*m_nodes[i].m_directionVector[2] - m_nodes[i].m_anchorDifference[2]*m_nodes[i].m_directionVector[0];
					if(i == 0 && writeTraces) std::cout << "scalarDirection: " << scalarDirection << "\n";
					if(i == 0 && writeTraces) std::cout << "scalarPerpendicular: " << scalarPerpendicular << "\n";
					m_nodes[i].m_collisionState = COLLIDE_COLLIDE_NOTANCHORCHANGE;
					if(scalarDirection > 0){
						frictionForceDirection.push_back(-springK[0][0]*scalarDirection*m_nodes[i].m_directionVector[0]);
						frictionForceDirection.push_back(0);
						frictionForceDirection.push_back(-springK[0][0]*scalarDirection*m_nodes[i].m_directionVector[2]);
						if(i == 0 && writeTraces) std::cout << "friction direction positive: " << frictionForceDirection[0] << " 0 " << frictionForceDirection[2] << std::endl;

						double frictionForceDirectionNorm = sqrt(frictionForceDirection[0]*frictionForceDirection[0]+frictionForceDirection[2]*frictionForceDirection[2]);

						if(i == 0 && writeTraces) std::cout << "check: " << frictionForceDirectionNorm << " " << mu[0][0]*m_nodes[i].m_penalty << std::endl;
						if(frictionForceDirectionNorm > mu[0][0]*m_nodes[i].m_penalty){
							m_nodes[i].m_collisionState = COLLIDE_COLLIDE_ANCHORCHANGE_1;
							frictionForceDirection[0] = mu[0][0]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[0];
							frictionForceDirection[2] = mu[0][0]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[2];
						}

						if(i == 0 && writeTraces) std::cout << "friction direction positive: " << frictionForceDirection[0] << " 0 " << frictionForceDirection[2] << std::endl;
						m_nodes[i].m_frictionElection = FORWARD;
						if(i == 0 && writeTraces) std::cout << "scalarDirection positive: " << mu[0][0]*scalarDirection << " " << scalarDirection << std::endl;
					}
					else{
						frictionForceDirection.push_back(-springK[0][0]*scalarDirection*m_nodes[i].m_directionVector[0]);
						frictionForceDirection.push_back(0);
						frictionForceDirection.push_back(-springK[0][0]*scalarDirection*m_nodes[i].m_directionVector[2]);
						if(i == 0 && writeTraces) std::cout << "friction direction negative: " << frictionForceDirection[0] << " 0 " << frictionForceDirection[2] << std::endl;

						double frictionForceDirectionNorm = sqrt(frictionForceDirection[0]*frictionForceDirection[0]+frictionForceDirection[2]*frictionForceDirection[2]);

						if(i == 0 && writeTraces) std::cout << "check: " << frictionForceDirectionNorm << " " << mu[0][1]*m_nodes[i].m_penalty << std::endl;
						if(frictionForceDirectionNorm > mu[0][1]*m_nodes[i].m_penalty){
							m_nodes[i].m_collisionState = COLLIDE_COLLIDE_ANCHORCHANGE_1;
							frictionForceDirection[0] = -mu[0][1]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[0];
							frictionForceDirection[2] = -mu[0][1]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[2];
						}
						if(i == 0 && writeTraces) std::cout << "friction direction negative: " << frictionForceDirection[0] << " 0 " << frictionForceDirection[2] << std::endl;
						m_nodes[i].m_frictionElection = BACKWARD;
						if(i == 0 && writeTraces) std::cout << "scalarDirection negative: " << mu[0][1]*scalarDirection << " " << scalarDirection << std::endl;
					}

					if(scalarPerpendicular > 0){
						frictionForcePerpendicular.push_back(-springK[0][1]*scalarPerpendicular*m_nodes[i].m_directionVector[2]);
						frictionForcePerpendicular.push_back(0);
						frictionForcePerpendicular.push_back(springK[0][1]*scalarPerpendicular*m_nodes[i].m_directionVector[0]);
						if(i == 0 && writeTraces) std::cout << "friction perpendicular positive: " << frictionForceDirection[0] << " 0 " << frictionForceDirection[2] << std::endl;
						
						double frictionForcePerpendicularNorm = sqrt(frictionForcePerpendicular[0]*frictionForcePerpendicular[0]+frictionForcePerpendicular[2]*frictionForcePerpendicular[2]);

						if(i == 0 && writeTraces) std::cout << "check: " << frictionForcePerpendicularNorm << " " << mu[0][2]*m_nodes[i].m_penalty << std::endl;
						if(frictionForcePerpendicularNorm > mu[0][2]*m_nodes[i].m_penalty){
							if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_NOTANCHORCHANGE){
								m_nodes[i].m_collisionState = COLLIDE_COLLIDE_ANCHORCHANGE_2;
							}
							else if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_ANCHORCHANGE_1){
								m_nodes[i].m_collisionState = COLLIDE_COLLIDE_ANCHORCHANGE_12;
							}
							frictionForcePerpendicular[0] = mu[0][2]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[2];
							frictionForcePerpendicular[2] = -mu[0][2]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[0];
						}

						if(i == 0 && writeTraces) std::cout << "friction perpendicular positive: " << frictionForceDirection[0] << " 0 " << frictionForceDirection[2] << std::endl;
						m_nodes[i].m_frictionElection2 = FORWARD;
						if(i == 0 && writeTraces) std::cout << "scalarPerpendicular positive: " << mu[0][0]*scalarPerpendicular << " " << scalarPerpendicular << std::endl;
					}
					else{
						frictionForcePerpendicular.push_back(-springK[0][1]*scalarPerpendicular*m_nodes[i].m_directionVector[2]);
						frictionForcePerpendicular.push_back(0);
						frictionForcePerpendicular.push_back(springK[0][1]*scalarPerpendicular*m_nodes[i].m_directionVector[0]);
						if(i == 0 && writeTraces) std::cout << "friction perpendicular negative: " << frictionForcePerpendicular[0] << " 0 " << frictionForcePerpendicular[2] << std::endl;

						double frictionForcePerpendicularNorm = sqrt(frictionForcePerpendicular[0]*frictionForcePerpendicular[0]+frictionForcePerpendicular[2]*frictionForcePerpendicular[2]);

						if(i == 0 && writeTraces) std::cout << "check: " << frictionForcePerpendicularNorm << " " << mu[0][2]*m_nodes[i].m_penalty << std::endl;
						if(frictionForcePerpendicularNorm > mu[0][2]*m_nodes[i].m_penalty){
							if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_NOTANCHORCHANGE){
								m_nodes[i].m_collisionState = COLLIDE_COLLIDE_ANCHORCHANGE_2;
							}
							else if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_ANCHORCHANGE_1){
								m_nodes[i].m_collisionState = COLLIDE_COLLIDE_ANCHORCHANGE_12;
							}
							frictionForcePerpendicular[0] = -mu[0][2]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[2];
							frictionForcePerpendicular[2] = mu[0][2]*m_nodes[i].m_penalty*m_nodes[i].m_directionVector[0];
						}
						if(i == 0 && writeTraces) std::cout << "friction perpendicular negative: " << frictionForcePerpendicular[0] << " 0 " << frictionForcePerpendicular[2] << std::endl;
						m_nodes[i].m_frictionElection2 = BACKWARD;
						if(i == 0 && writeTraces) std::cout << "scalarPerpendicular negative: " << mu[0][1]*scalarPerpendicular << " " << scalarDirection << std::endl;
					}

					
					auxForce[0] = frictionForceDirection[0] + frictionForcePerpendicular[0];
					auxForce[2] = frictionForceDirection[2] + frictionForcePerpendicular[2];

				
					if(m_nodes[i].m_collisionState != COLLIDE_COLLIDE_NOTANCHORCHANGE){
						//update anchor
						Eigen::Matrix2f A;
						Eigen::Vector2f b;

						if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_ANCHORCHANGE_1){
							A(0,0)=-m_nodes[i].m_directionVector[0];
							A(0,1)=-m_nodes[i].m_directionVector[2];
							A(1,0)=-m_nodes[i].m_directionVector[2];
							A(1,1)= m_nodes[i].m_directionVector[0];
							if(m_nodes[i].m_frictionElection == FORWARD){
								b(0) = -mu[0][0]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = scalarPerpendicular - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}
							else{
								b(0) = mu[0][1]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = scalarPerpendicular - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}

						}
						else if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_ANCHORCHANGE_2){
							A(0,0)=-m_nodes[i].m_directionVector[0];
							A(0,1)=-m_nodes[i].m_directionVector[2];
							A(1,0)=-m_nodes[i].m_directionVector[2];
							A(1,1)= m_nodes[i].m_directionVector[0];
							if(m_nodes[i].m_frictionElection2 == FORWARD){
								b(0) = scalarDirection - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = -mu[0][2]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}
							else{
								b(0) = scalarDirection - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = mu[0][2]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}
						}
						else if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_ANCHORCHANGE_12){
							A(0,0)=-m_nodes[i].m_directionVector[0];
							A(0,1)=-m_nodes[i].m_directionVector[2];
							A(1,0)=-m_nodes[i].m_directionVector[2];
							A(1,1)= m_nodes[i].m_directionVector[0];
							if(m_nodes[i].m_frictionElection == FORWARD && m_nodes[i].m_frictionElection2 == FORWARD){
								b(0) = -mu[0][0]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = -mu[0][2]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}
							else if(m_nodes[i].m_frictionElection == FORWARD && m_nodes[i].m_frictionElection2 == BACKWARD){
								b(0) = -mu[0][0]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = mu[0][2]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}
							else if(m_nodes[i].m_frictionElection == BACKWARD && m_nodes[i].m_frictionElection2 == FORWARD){
								b(0) = mu[0][1]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = -mu[0][2]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}
							else if(m_nodes[i].m_frictionElection == BACKWARD && m_nodes[i].m_frictionElection2 == BACKWARD){
								b(0) = mu[0][1]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[0] - x[i][0]*m_nodes[i].m_directionVector[2];
								b(1) = mu[0][2]*m_nodes[i].m_penalty/springK[0][0] - x[i][0]*m_nodes[i].m_directionVector[2] + x[i][0]*m_nodes[i].m_directionVector[0];
							}
						}

						Eigen::Vector2f res = A.fullPivLu().solve(b);

						m_nodes[i].m_anchorDifference[0]=res(0);
						m_nodes[i].m_anchorDifference[2]=res(2);

						//check

					}
				}

				f[i][0] += auxForce[0];
				f[i][1] += auxForce[1];
				f[i][2] += auxForce[2];

				frictionForceSum[0] += f[i][0];
				frictionForceSum[1] += f[i][1];
				frictionForceSum[2] += f[i][2];

				//if(i == 28 & writeTraces) std::cout << "friction: " << getContext()->getTime() << " " << f[i] << std::endl;
			}
		}
		else{
			if(m_nodes[i].m_collisionState < 5){
				m_nodes[i].m_collisionState = COLLIDE_NOTCOLLIDE;
			}
			else{
				m_nodes[i].m_collisionState = NOTCOLLIDE_NOTCOLLIDE;
			}
		}
	}

	if(writeTraces) std::cout << "####sum: " << frictionForceSum[0] << " " << frictionForceSum[1] << " " << frictionForceSum[2] << "#######" << std::endl;
}

template <class DataTypes> 
void Friction<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
	//init
	VecDeriv& f = *d_df.beginEdit();
	const VecDeriv& dx = d_dx.getValue();
	const VecCoord& mu = m_mu.getValue();
	const VecCoord& springK = m_springK.getValue();
	bool writeTraces = f_writeTracesD.getValue();

	//std::vector<sofa::component::linearsolver::MyCGLinearSolver< sofa::component::linearsolver::GraphScatteredMatrix, sofa::component::linearsolver::GraphScatteredVector> * > sf;
	//this->getContext()->get<sofa::component::linearsolver::MyCGLinearSolver< sofa::component::linearsolver::GraphScatteredMatrix, sofa::component::linearsolver::GraphScatteredVector> >(&sf,core::objectmodel::BaseContext::SearchUp);
				
	//I_3
	Eigen::MatrixXd I3 = Eigen::MatrixXd::Identity(3,3);

	//I_3 - n*n^t
	Eigen::MatrixXd I3_NNt(3, 3);
	I3_NNt(0,0) = 1;I3_NNt(0,1) = 0;I3_NNt(0,2) = 0;I3_NNt(1,0) = 0;I3_NNt(1,1) = 0;I3_NNt(1,2) = 0;
	I3_NNt(2,0) = 0;I3_NNt(2,1) = 0;I3_NNt(2,2) = 1;

	//n*n^t
	Eigen::MatrixXd NNt(3, 3);
	NNt(0,0) = 0;NNt(0,1) = 0;NNt(0,2) = 0;NNt(1,0) = 0;NNt(1,1) = 1;NNt(1,2) = 0;
	NNt(2,0) = 0;NNt(2,1) = 0;NNt(2,2) = 0;
	//everything simplified to normal vector n=(0,1,0)

	for(unsigned int i = 0; i < f.size(); i++){
		if(!m_frictionRelations[i]) continue;
		

		if(i == 0 && writeTraces)std::cout << "initial f: " << f[i] << std::endl;
		//penalty for all of them
		f[i][1] += -m_floorK.getValue()*dx[i][1]*mparams->kFactor();
		/*if(i == 0 && writeTraces)std::cout << "after penalty: " << f[i] << std::endl;
		Eigen::MatrixXd auxp(3, 1);
		auxp(0,0) = dx[i][0];auxp(1,0) = dx[i][1];auxp(2,0) = dx[i][2];
		if(i == 0 && writeTraces)std::cout << "auxp: " << auxp << std::endl;

		if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_NOTANCHORCHANGE){
			Eigen::MatrixXd result(3, 3);
			result = I3_NNt * (-springK[0][0]);

			Eigen::MatrixXd aux(3, 1);
			aux = result * auxp;

			f[i][0] += aux(0,0)*mparams->kFactor();
			f[i][1] += aux(1,0)*mparams->kFactor();
			f[i][2] += aux(2,0)*mparams->kFactor();

			if(i == 0 && writeTraces)std::cout << "result: " << result << std::endl;
			if(i == 0 && writeTraces)std::cout << "aux: " << aux << std::endl;
			if(i == 0 && writeTraces)std::cout << "k factor: " << mparams->kFactor() << std::endl;
			if(i == 0 && writeTraces)std::cout << "final f: " << f[i] << std::endl;
		}

		else if(m_nodes[i].m_collisionState == COLLIDE_COLLIDE_ANCHORCHANGE_1){

			Eigen::MatrixXd auxq(3, 1);
			auxq(0,0) = dx[m_frictionRelations[i]][0];auxq(1,0) = dx[m_frictionRelations[i]][1];auxq(2,0) = dx[m_frictionRelations[i]][2];
			if(i == 0 && writeTraces)std::cout << "auxq: " << auxq << std::endl;

			//F_n^t
			Eigen::MatrixXd Fnt(1,3);
			Fnt(0,0) = 0;Fnt(0,1) = m_nodes[i].m_penalty;Fnt(0,2) = 0;

			//p-a
			Eigen::MatrixXd P_Au(3,1);
			P_Au(0,0) = m_nodes[i].m_anchorDifferenceUnit[0];P_Au(1,0) = 0;P_Au(2,0) = m_nodes[i].m_anchorDifferenceUnit[2];

			//(p - a)^t
			Eigen::MatrixXd P_Aut(1,3);
			P_Aut = P_Au.transpose();

			//d(p-a)/dp (unit)
			Eigen::MatrixXd Dp_aDp(3,3);
			Dp_aDp = (I3 - P_Au*P_Aut) * I3_NNt * (1/m_nodes[i].m_differenceNorm);

			//q - p
			Eigen::MatrixXd q_p(3,1);
			q_p(0,0)=m_nodes[i].m_directionVector[0]; q_p(1,0)=0; q_p(2,0)=m_nodes[i].m_directionVector[2];

			//q - p transposed
			Eigen::MatrixXd q_pt(1,3);
			q_pt(0,0)=m_nodes[i].m_directionVector[0]; q_pt(0,1)=0; q_pt(0,2)=m_nodes[i].m_directionVector[2];

			//q - p perpendicular
			Eigen::MatrixXd q_p_p(3,1);
			q_p_p(0,0)=m_nodes[i].m_directionVector[2]; q_p_p(1,0)=0; q_p_p(2,0)=-m_nodes[i].m_directionVector[0];

			//q - p perpendicular transposed
			Eigen::MatrixXd q_p_pt(1,3);
			q_p_pt(0,0)=m_nodes[i].m_directionVector[2]; q_p_pt(0,1)=0; q_p_pt(0,2)=-m_nodes[i].m_directionVector[0];

			//q - p differentiated
			Eigen::MatrixXd Dq_pDp(3,3);
			Dq_pDp = (I3 - q_p*q_pt) * I3_NNt * (-1/m_nodes[i].m_directionNorm);

			//q - p differentiated (in q)
			Eigen::MatrixXd Dq_pDq(3,3);
			Dq_pDq = (I3 - q_p*q_pt) * I3_NNt * (1/m_nodes[i].m_directionNorm);

			//q - p perpendicular differentiated
			Eigen::MatrixXd Dq_p_pDp(3,3);
			Dq_p_pDp = (I3 - q_p_p*q_p_pt) * I3_NNt * (-1/m_nodes[i].m_directionNorm);

			//q - p perpendicular differentiated (in q)
			Eigen::MatrixXd Dq_p_pDq(3,3);
			Dq_p_pDq = (I3 - q_p_p*q_p_pt) * I3_NNt * (1/m_nodes[i].m_directionNorm);

			Eigen::MatrixXd result(3, 3);
			Eigen::MatrixXd result2(3, 3);

			if(m_nodes[i].m_frictionElection == FORWARD && m_nodes[i].m_frictionElection == FORWARD){
				Eigen::MatrixXd Dmu_Dp(1, 3);
				Dmu_Dp = (((P_Aut * Dq_pDp + q_pt*Dp_aDp)*(mu[0][0]) + (P_Aut * Dq_p_pDp + q_p_pt*Dp_aDp)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (P_Aut * Dq_pDp + q_pt*Dp_aDp + P_Aut * Dq_p_pDp + q_p_pt*Dp_aDp) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result = Dp_aDp * I3_NNt * (m_nodes[i].m_currentMu*m_nodes[i].m_penalty) + 
					P_Au * Fnt * NNt (m_nodes[i].m_currentMu*m_floorK.getValue()/m_nodes[i].m_penalty) +
					P_Au * Dmu_Dp * (m_nodes[i].m_penalty);

				Eigen::MatrixXd Dmu_Dq(1, 3);
				Dmu_Dq = (((P_Aut * Dq_pDq)*(mu[0][0]) + (P_Aut * Dq_p_pDq)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (P_Aut * Dq_pDq + P_Aut * Dq_p_pDq) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result2 = P_Au * Dmu_Dq * (m_nodes[i].m_penalty);
			}
			else if(m_nodes[i].m_frictionElection == FORWARD && m_nodes[i].m_frictionElection == BACKWARD){
				Eigen::MatrixXd Dmu_Dp(1, 3);
				Dmu_Dp = (((P_Aut * Dq_pDp + q_pt*Dp_aDp)*(mu[0][0]) - (P_Aut * Dq_p_pDp + q_p_pt*Dp_aDp)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (P_Aut * Dq_pDp + q_pt*Dp_aDp - P_Aut * Dq_p_pDp - q_p_pt*Dp_aDp) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result = Dp_aDp * I3_NNt * (m_nodes[i].m_currentMu*m_nodes[i].m_penalty) + 
					P_Au * Fnt * NNt (m_nodes[i].m_currentMu*m_floorK.getValue()/m_nodes[i].m_penalty) +
					P_Au * Dmu_Dp * (m_nodes[i].m_penalty);

				Eigen::MatrixXd Dmu_Dq(1, 3);
				Dmu_Dq = (((P_Aut * Dq_pDq)*(mu[0][0]) - (P_Aut * Dq_p_pDq)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (P_Aut * Dq_pDq - P_Aut * Dq_p_pDq) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result2 = P_Au * Dmu_Dq * (m_nodes[i].m_penalty);
			}
			else if(m_nodes[i].m_frictionElection == BACKWARD && m_nodes[i].m_frictionElection == FORWARD){
				Eigen::MatrixXd Dmu_Dp(1, 3);
				Dmu_Dp = (((-P_Aut * Dq_pDp - q_pt*Dp_aDp)*(mu[0][1]) + (P_Aut * Dq_p_pDp + q_p_pt*Dp_aDp)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (-P_Aut * Dq_pDp - q_pt*Dp_aDp + P_Aut * Dq_p_pDp + q_p_pt*Dp_aDp) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result = Dp_aDp * I3_NNt * (m_nodes[i].m_currentMu*m_nodes[i].m_penalty) + 
					P_Au * Fnt * NNt (m_nodes[i].m_currentMu*m_floorK.getValue()/m_nodes[i].m_penalty) +
					P_Au * Dmu_Dp * (m_nodes[i].m_penalty);

				Eigen::MatrixXd Dmu_Dq(1, 3);
				Dmu_Dq = (((-P_Aut * Dq_pDq)*(mu[0][1]) + (P_Aut * Dq_p_pDq)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (-P_Aut * Dq_pDq + P_Aut * Dq_p_pDq) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result2 = P_Au * Dmu_Dq * (m_nodes[i].m_penalty);
			}
			else if(m_nodes[i].m_frictionElection == BACKWARD && m_nodes[i].m_frictionElection == BACKWARD){
				Eigen::MatrixXd Dmu_Dp(1, 3);
				Dmu_Dp = (((-P_Aut * Dq_pDp - q_pt*Dp_aDp)*(mu[0][1]) - (P_Aut * Dq_p_pDp + q_p_pt*Dp_aDp)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (-P_Aut * Dq_pDp - q_pt*Dp_aDp - P_Aut * Dq_p_pDp - q_p_pt*Dp_aDp) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result = Dp_aDp * I3_NNt * (m_nodes[i].m_currentMu*m_nodes[i].m_penalty) + 
					P_Au * Fnt * NNt (m_nodes[i].m_currentMu*m_floorK.getValue()/m_nodes[i].m_penalty) +
					P_Au * Dmu_Dp * (m_nodes[i].m_penalty);

				Eigen::MatrixXd Dmu_Dq(1, 3);
				Dmu_Dq = (((-P_Aut * Dq_pDq)*(mu[0][1]) - (P_Aut * Dq_p_pDq)*(mu[0][2])) *  m_nodes[i].m_currentScalarSum - (-P_Aut * Dq_pDq - P_Aut * Dq_p_pDq) * (m_nodes[i].m_currentScalarMuSum)) * (1/m_nodes[i].m_currentScalarSum*1/m_nodes[i].m_currentScalarSum);
				result2 = P_Au * Dmu_Dq * (m_nodes[i].m_penalty);
			}

			Eigen::MatrixXd aux(3, 1);
			aux = result * auxp;

			f[i][0] += aux(0,0)*mparams->kFactor();
			f[i][1] += aux(1,0)*mparams->kFactor();
			f[i][2] += aux(2,0)*mparams->kFactor();

			if(i == 0 && writeTraces)std::cout << "result: " << result << std::endl;
			if(i == 0 && writeTraces)std::cout << "aux: " << aux << std::endl;
			if(i == 0 && writeTraces)std::cout << "k factor: " << mparams->kFactor() << std::endl;
			if(i == 0 && writeTraces)std::cout << "final f: " << f[i] << std::endl;

			Eigen::MatrixXd aux2(3, 1);
			aux2 = result2 * auxq;

			if(i == 0 && writeTraces)std::cout << "initial f_q: " << f[m_frictionRelations[i]] << std::endl;

			f[m_frictionRelations[i]][0] += aux2(0,0)*mparams->kFactor();
			f[m_frictionRelations[i]][1] += aux2(1,0)*mparams->kFactor();
			f[m_frictionRelations[i]][2] += aux2(2,0)*mparams->kFactor();

			if(i == 0 && writeTraces)std::cout << "result2: " << result2 << std::endl;
			if(i == 0 && writeTraces)std::cout << "aux2: " << aux2 << std::endl;
			if(i == 0 && writeTraces)std::cout << "k factor: " << mparams->kFactor() << std::endl;
			if(i == 0 && writeTraces)std::cout << "final f_q: " << f[m_frictionRelations[i]] << std::endl;*/

	}
}

template<class DataTypes>
void Friction<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
	std::cout << "*********************************************************************PPPPPPPPPPPPPPP" << std::endl;

}


		} // namespace forcefield
	} // namespace component
} // namespace sofa

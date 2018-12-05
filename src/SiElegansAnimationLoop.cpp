#define SOFA_COMPONENT_FORCEFIELD_SiElegansAnimationLoop_CPP

#include "SiElegansAnimationLoop.h"


#include <sofa/core/ObjectFactory.h>

#include <sofa/simulation/common/PrintVisitor.h>
#include <sofa/simulation/common/FindByTypeVisitor.h>
#include <sofa/simulation/common/ExportGnuplotVisitor.h>
#include <sofa/simulation/common/InitVisitor.h>
#include <sofa/simulation/common/AnimateVisitor.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/simulation/common/CollisionVisitor.h>
//#include <sofa/simulation/common/CollisionBeginEvent.h>
//#include <sofa/simulation/common/CollisionEndEvent.h>
#include <sofa/simulation/common/UpdateContextVisitor.h>
#include <sofa/simulation/common/UpdateMappingVisitor.h>
#include <sofa/simulation/common/ResetVisitor.h>
#include <sofa/simulation/common/VisualVisitor.h>
#include <sofa/simulation/common/ExportOBJVisitor.h>
#include <sofa/simulation/common/WriteStateVisitor.h>
#include <sofa/simulation/common/XMLPrintVisitor.h>
#include <sofa/simulation/common/PropagateEventVisitor.h>
#include <sofa/simulation/common/BehaviorUpdatePositionVisitor.h>
#include <sofa/simulation/common/AnimateBeginEvent.h>
#include <sofa/simulation/common/AnimateEndEvent.h>
#include <sofa/simulation/common/UpdateMappingEndEvent.h>
#include <sofa/simulation/common/CleanupVisitor.h>
#include <sofa/simulation/common/DeleteVisitor.h>
#include <sofa/simulation/common/UpdateBoundingBoxVisitor.h>
#include <sofa/simulation/common/xml/NodeElement.h>

#include <sofa/helper/system/SetDirectory.h>
//#include <sofa/helper/system/PipeProcess.h>
#include <sofa/helper/AdvancedTimer.h>

#include <sofa/core/visual/VisualParams.h>

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <curl/curl.h>
#include <curl/easy.h>
#include <string>
#include <sstream>
#include <vector>




std::string buffer;

namespace sofa
{

namespace simulation
{


SOFA_DECL_CLASS(SiElegansAnimationLoop);

int SiElegansAnimationLoopClass = core::RegisterObject("Si Elegans animation loop")
        .add< SiElegansAnimationLoop >()
        ;



SiElegansAnimationLoop::SiElegansAnimationLoop(simulation::Node* _gnode)
    : Inherit()
  , f_getURL(initData(&f_getURL,std::string("http://127.0.0.1:3000/PE/api/input"),"URL_get","URL for getting muscle output"))
  , f_finishURL(initData(&f_finishURL,std::string("http://10.0.20.26:1730/PEIF/api/simulation_end/sofa"),"URL_finish","URL for finishing the simulation"))
  , f_postURL(initData(&f_postURL,std::string("http://127.0.0.1:3000/PE/api/runtime_stimuli"),"URL_post","URL for posting sensory input"))
  , f_writeTraces( initData(&f_writeTraces, true, "write_traces","write traces in console") )
  , f_finishingTime(initData(&f_finishingTime, 0, "finishing_time", "time when the finishing curl get will be done") )
  , f_syncTime(initData(&f_syncTime, 20, "sync_time", "how often PE syncs with neurons") )
{
    //assert(gnode);
  UNUSED(_gnode);
}

SiElegansAnimationLoop::~SiElegansAnimationLoop()
{

}

void SiElegansAnimationLoop::init()
{
  std::cout << "init Loop" << std::endl;
  Inherit::init();
  sofa::helper::vector<component::behaviormodel::Stimuli*> aux;
  gnode->getTreeObjects<component::behaviormodel::Stimuli>(&aux);
  m_stimuliVec = aux;
  if (f_writeTraces.getValue()){
    for (unsigned int i=0; i<m_stimuliVec.size(); i++)
	std::cout << "added stimulus: " << m_stimuliVec[i]->getName() << std::endl;
  }
  
  sofa::helper::vector<sofa::component::mass::InterfaceManagerActivation<defaulttype::Vec3dTypes,double>*> aux2;
  gnode->getTreeObjects<sofa::component::mass::InterfaceManagerActivation<defaulttype::Vec3dTypes,double>>(&aux2);
  if(aux2.size() > 0){
	m_muscleActivation = aux2.at(0);
  }
  else{
	  m_muscleActivation = NULL;
  }
  m_currentSyncTime = float(float(f_syncTime.getValue())/1000);
}

std::vector<std::string> &SiElegansAnimationLoop::split(std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> SiElegansAnimationLoop::split(std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

size_t writefunc(void *ptr, size_t size, size_t nmemb, std::string *s)
{
  UNUSED(s);
  buffer = (char*)ptr;  
  return size*nmemb;
}

void SiElegansAnimationLoop::step(const core::ExecParams* params, SReal dt)
{
  bool writeTraces = f_writeTraces.getValue();
  if(writeTraces) std::cout <<  "************step************  " << gnode->getTime() << std::endl;
  
  SReal startTime = gnode->getTime();
  if (dt == 0) dt = this->gnode->getDt();
  
  CURL *curl;
  CURLcode res;

  bool inputRead = false;
  
  // In windows, this will init the winsock stuff 
  curl_global_init(CURL_GLOBAL_ALL);
  
  curl = curl_easy_init();
  std::string getURL = f_getURL.getValue();
  std::vector<double> muscleValues;
 //std::cout << "coño: " << gnode->getTime() << " " <<  float(float(f_finishingTime.getValue())/1000)-0.001 << std::endl;
  if(gnode->getTime() >= float(float(f_finishingTime.getValue())/1000)-0.001){
	if(writeTraces) std::cout << "finishing GET" << std::endl;
	  if(curl && f_finishURL.getValue() != "no_CURL") {
		  curl_easy_setopt(curl, CURLOPT_URL, f_finishURL.getValue().c_str());
		  curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
		  curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writefunc);
		  res = curl_easy_perform(curl);
		return;
	  }
  }

  bool syncNeurons = false;
  if(writeTraces) std::cout << "time: " << m_currentSyncTime << " " << startTime << "\n";
  if(m_currentSyncTime > float(float(f_syncTime.getValue())/1000)-0.001){
	  m_currentSyncTime -= float(f_syncTime.getValue())/1000;
	  m_currentSyncTime += dt;
	  syncNeurons = true;
	  if(writeTraces) std::cout << "sync " << m_currentSyncTime << "\n";
  }
  else{
	  m_currentSyncTime += dt;
  }

  if(curl && getURL != "no_CURL" && syncNeurons) {

	  if(writeTraces) std::cout << "get\n";
    while(!inputRead){
	//std::cout << "input GET" << std::endl;
      curl_easy_setopt(curl, CURLOPT_URL, getURL.c_str());
      curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writefunc);
      res = curl_easy_perform(curl);
      
      if(res!=CURLE_OK){
        fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
      }
      else{
        std::string the_output = buffer.c_str();
        if(the_output.substr(0,15) == "No muscles data"){
          //if(writeTraces) std::cout <<  "No muscles data \n" ;          
        }
        else{	
		std::vector<std::string> values = split(the_output, ',');
         	muscleValues.clear();
          	for(unsigned int i = 0; i < values.size(); i++){
            		//if(writeTraces) std::cout << i << ": " << values.at(i).c_str() << std::endl;
            		//muscleValues.push_back(stod(values.at(i)));
			double d = atof(values.at(i).c_str());
			double ds = stod(values.at(i));
			// std::cout <<"----------- "<< d << ": " << ds << std::endl;			
            		muscleValues.push_back(d);
	    	}
		if(m_muscleActivation){
			m_muscleActivation->setMuscleActivation(muscleValues);
			inputRead = true;
		}
        }
      }
    }
    // always cleanup  
    curl_easy_cleanup(curl);
  }

    sofa::helper::AdvancedTimer::stepBegin("AnimationStep");

    sofa::helper::AdvancedTimer::begin("Animate");

#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printNode("Step");
#endif

    {
        AnimateBeginEvent ev ( dt );
        PropagateEventVisitor act ( params, &ev );
        gnode->execute ( act );
    }

    
    BehaviorUpdatePositionVisitor beh(params , dt);
    gnode->execute ( beh );
    AnimateVisitor act(params, dt);
    gnode->execute ( act );

    gnode->setTime ( startTime + dt );
    gnode->execute< UpdateSimulationContextVisitor >(params);

    {
        AnimateEndEvent ev ( dt );
        PropagateEventVisitor act ( params, &ev );
        gnode->execute ( act );
    }

    sofa::helper::AdvancedTimer::stepBegin("UpdateMapping");
    //Visual Information update: Ray Pick add a MechanicalMapping used as VisualMapping
    gnode->execute< UpdateMappingVisitor >(params);
    sofa::helper::AdvancedTimer::step("UpdateMappingEndEvent");
    {
        UpdateMappingEndEvent ev ( dt );
        PropagateEventVisitor act ( params , &ev );
        gnode->execute ( act );
    }
    sofa::helper::AdvancedTimer::stepEnd("UpdateMapping");
    

#ifndef SOFA_NO_UPDATE_BBOX
    sofa::helper::AdvancedTimer::stepBegin("UpdateBBox");
    gnode->execute< UpdateBoundingBoxVisitor >(params);
    sofa::helper::AdvancedTimer::stepEnd("UpdateBBox");
#endif
#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printCloseNode("Step");
#endif

    for (std::vector<component::behaviormodel::Stimuli*>::iterator it=m_stimuliVec.begin(); it != m_stimuliVec.end(); ++it)
      (*it)->step();
   
    std::string currentStimuli;
    for (std::vector<component::behaviormodel::Stimuli*>::iterator it=m_stimuliVec.begin(); it != m_stimuliVec.end(); ++it)
    {
      component::behaviormodel::Stimuli* stimulus = *it;
      stimulus->createStimuliString();
      currentStimuli += stimulus->getStimuliString();
      stimulus->clearStimuliString();
      stimulus->resetStimuli();
    }
    //if (writeTraces) std::cout <<  "--- stimuli output str:" << currentStimuli << std::endl;


	//code to write the NaCl stimuli in ASEL and ASER in a file
	/*std::ofstream out("C:/Users/amujika/Box Sync/Kodigo/sielegans-bitbucket/PhysicsEngine/physics-engine/sofa/applications/plugins/SiElegansPlugin/Examples/data/alicia.txt", std::ios::app);
	//write stimuli for Alicia. ASEL = 40, ASER = 41, NaCl = 65535
	std::vector<std::string> stms = split(currentStimuli, '#');
	for(unsigned int i = 0; i < stms.size(); i++){
		std::vector<std::string> stm = split(stms.at(i), ',');
		if(stm.size() == 1){
			out << "t: " << stm.at(0) << '\n';
		}
		else if(stm.size() == 3){
			if(stm.at(0) == "40" && stm.at(1)=="65535" ){
				out << "ASEL, NaCl: " << stm.at(2) << "\n";
			}
			if(stm.at(0) == "41" && stm.at(1)=="65535" ){
				out << "ASER, NaCl: " << stm.at(2) << "\n";
			}
		}
	}*/

   
	//POST
    curl = curl_easy_init();
    std::string postURL = f_postURL.getValue();

    if(curl && postURL != "no_CURL"  && syncNeurons) {
		std::cout << "post\n";
		if(writeTraces) std::cout << "stimuli POST" << std::endl;
      struct curl_slist *headers = NULL;
      headers = curl_slist_append(headers, "Content-Type: application/json");
  
      curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
    curl_easy_setopt(curl, CURLOPT_USERNAME,"admin");
    curl_easy_setopt(curl, CURLOPT_PASSWORD, "1sielegans!");
      curl_easy_setopt(curl, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
      // dirección real del peif: ip_server/PE/api/input
      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, false);
      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, false);
      curl_easy_setopt(curl, CURLOPT_URL, postURL.c_str());
    //curl_easy_setopt(curl, CURLOPT_URL, "https://150.241.250.4:2730/PE/api/input");
      // example.com is redirected, so we tell libcurl to follow redirection 
      curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
      // Perform the request, res will get the return code 
      curl_easy_setopt(curl, CURLOPT_VERBOSE, 1);
      // example.com is redirected, so we tell libcurl to follow redirection
      //curl_easy_setopt(curl, CURLOPT_UPLOAD, 1L);
      curl_easy_setopt(curl, CURLOPT_CUSTOMREQUEST, "POST");
      //curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);
      curl_easy_setopt(curl, CURLOPT_POSTFIELDS, currentStimuli.c_str());
      // Perform the request, res will get the return code
      res = curl_easy_perform(curl);
      //std::cout << res << std::endl;
      // Check for errors 
      if(res != CURLE_OK)
  fprintf(stderr, "curl_easy_perform() failed: %s\n",
    curl_easy_strerror(res));
  
      // always cleanup 
      curl_easy_cleanup(curl);
    }
   
    ///////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////

    sofa::helper::AdvancedTimer::end("Animate");
    sofa::helper::AdvancedTimer::stepEnd("AnimationStep");
}


} // namespace simulation

} // namespace sofa

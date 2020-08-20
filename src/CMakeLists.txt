include(${SOFA_CMAKE_DIR}/preProject.cmake)
project(SiElegans VERSION 1.0)

find_package(CURL)
find_package(Qt4)

include_directories(
	${CURL_INCLUDE_DIR}
	${QT_INCLUDE_DIR}
)

set(HEADER_FILES

    SiElegansAnimationLoop.h
	MuscleActivation.h
	MuscleActivation.inl
	initPlugin.h
	StandardTetrahedralFEMForceFieldTwoMaterials.inl
	StandardTetrahedralFEMForceFieldTwoMaterials.h
	Friction.h
	Friction.inl
	CGLinearSolverTraces.h
	CGLinearSolverTraces.inl
	Stimuli.h StimuliChemo.h StimuliMechano.h  StimuliProprio.h StimuliThermo.h
	StimuliBarycentricPenalityContact.h
	StimuliBarycentricPenalityContact.inl
	StimuliPenalityContactForceField.h
	StimuliPenalityContactForceField.inl
	PatternPlusFileActivation.h
	PatternPlusFileActivation.inl
	

	)

set(SOURCE_FILES

    SiElegansAnimationLoop.cpp
	MuscleActivation.cpp
    initPlugin.cpp
	StandardTetrahedralFEMForceFieldTwoMaterials.cpp
	Friction.cpp
	CGLinearSolverTraces.cpp
	Stimuli.cpp StimuliChemo.cpp StimuliMechano.cpp  StimuliProprio.cpp StimuliThermo.cpp
	StimuliBarycentricPenalityContact.cpp
	StimuliPenalityContactForceField.cpp
	PatternPlusFileActivation.cpp

	)

set(README_FILES

	SiElegans.txt
	)

	SET( GUI_MOC_HDRS SaveScreenShots.h)

# Generate moc files
SOFA_QT4_WRAP_CPP(MOC_FILES ${GUI_MOC_HDRS})
#SOFA_QT4_WRAP_UI(FORM_FILES ${UI_FILES})


	
ToRelativePath(MOC_FILES "${CMAKE_CURRENT_SOURCE_DIR}" "${MOC_FILES}")
set(MOC_FILES_GROUP "Generated")
	
add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${MOC_FILES} ${SOURCE_FILES} ${README_FILES})

include_directories("include/rapidjson")
AddCompilerDefinitions("SOFA_BUILD_SIELEGANSPLUGIN")
AddSourceDependencies(eigen ${CURL_INCLUDE_DIR} ${QT_INCLUDE_DIR})
AddLinkerDependencies(${CURL_LIBRARY} SofaMiscForceField SofaDeformable SofaBoundaryCondition SofaMiscTopology SofaBaseMechanics SofaBaseTopology SofaBaseLinearSolver SofaSimulationTree SofaSimulationCommon SofaCore SofaDefaultType SofaHelper SofaRigid SofaOpenglVisual newmat SofaBaseVisual tinyxml SofaComponentCommon SofaComponentBase SofaComponentAdvanced SofaComponentGeneral SofaComponentMain SofaComponentMisc ${SOFA_VISUAL_LIB} SofaGuiQt)

include(${SOFA_CMAKE_DIR}/postProject.cmake)

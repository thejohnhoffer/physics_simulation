


#include "StandardTetrahedralFEMForceFieldTwoMaterials.h"

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
StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::StandardTetrahedralFEMForceFieldTwoMaterials() 
	: f_secondParameterSet(initData(&f_secondParameterSet,"SecondParameterSet","The global parameters specifying the material"))
	, f_secondMaterialStart( initData(&f_secondMaterialStart,(unsigned)0,"secondMaterialStart","secondMaterialStart") )
	, f_writeTraces(initData(&f_writeTraces, false, "write_traces", "write traces to debug" ))
{
}

template <class DataTypes> 
StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::~StandardTetrahedralFEMForceFieldTwoMaterials()
{
}


template<class DataType>
string StandardTetrahedralFEMForceFieldTwoMaterials<DataType>::intToString(int number)
{
	ostringstream convert;
	string numberSTR;
	
	convert << number;
	numberSTR = convert.str();
	convert.str("");

	return numberSTR;
}

template<class DataType>
int StandardTetrahedralFEMForceFieldTwoMaterials<DataType>::loadVector(string filename, int i)
{
	string line;
	int numLine = 0;
	vector<double> aux;

	ifstream file(filename);

	if(file.is_open())
	{
		while(getline(file, line))
		{
			numLine++;

			if(numLine == 1)
				m = std::stoi(line); // Rows

			else if(numLine == 2)
				n = std::stoi(line); // Columns

			else
				aux.push_back(std::stod(line));
		}
	}
	else
	{
		cout << "\nCan't open file\n" << endl;
		return -2;
	}

	file.close();

	internalForce.push_back(aux);

	return 0;
}

template<class DataType>
vector<double> StandardTetrahedralFEMForceFieldTwoMaterials<DataType>::interpolate(int row)
{
	vector<double> internalForceInterpolate;
	double f1, f2, fm;
	
	double deltaT = PERIOD/8.0f;
	double time = row * PERIOD / 8.0f;

	for(int i=0; i<internalForce[row].size(); i++)
	{
		f1 = internalForce[row][i];
		f2 = internalForce[row+1][i];

		fm = f1 + ((f2-f1)/(deltaT)) * (timestep - time);

		internalForceInterpolate.push_back(fm);
	}

	return internalForceInterpolate;
}


template<class DataTypes>
void StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    //const VecCoord& x = *this->mstate->getX();
    /*const VecCoord& x = this->mstate->readPositions().ref();
    unsigned int i=0,j=0,k=0,l=0;
    unsigned int nbEdges=this->_topology->getNbEdges();
    const vector< Edge> &edgeArray=this->_topology->getEdges() ;

    helper::vector<EdgeInformation>& edgeInf = *(this->edgeInfo.beginEdit());
    tetrahedronRestInfoVector& tetrahedronInf = *(this->tetrahedronInfo.beginEdit());

    EdgeInformation *einfo;

    TetrahedronRestInformation *tetInfo;
    unsigned int nbTetrahedra=this->_topology->getNbTetrahedra();
    const std::vector< Tetrahedron> &tetrahedronArray=this->_topology->getTetrahedra() ;

    for(i=0; i<nbTetrahedra; i+=1 )
    {
        tetInfo=&tetrahedronInf[i];
//            Matrix3 &df=tetInfo->deformationGradient;
//            Matrix3 Tdf=df.transposed();
        BaseMeshTopology::EdgesInTetrahedron te=this->_topology->getEdgesInTetrahedron(i);

        /// describe the jth vertex index of triangle no i 
        const Tetrahedron &ta= tetrahedronArray[i];
        for(j=0;j<6;j++)
		{
            einfo= &edgeInf[te[j]];
            Edge e=this->_topology->getLocalEdgesInTetrahedron(j);

            k=e[0];
            l=e[1];


                
            std::vector< Vector3 >  points;
            Coord pa = x[ta[k]];
            Coord pb = x[ta[l]];
            points.push_back(pa);
            points.push_back(pb);
            Vec<4,float> color2 = Vec<4,float>(1.0f , 0.0f , 0.0f,1.0f);
            vparams->drawTool()->drawLines(points,1,color2 );

        }// end of for j                
    }//end of for i


    

    if (!vparams->displayFlags().getShowForceFields()) return;
    if (!this->mstate) return;

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    if (vparams->displayFlags().getShowWireFrame())
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);*/

	int count = 0;
	const VecCoord& x = this->mstate->readPositions().ref();

	if(m_frictionRelations->size() > 0){
		for(unsigned int i = 0; i < m_nodes->size(); i++){
			if(m_frictionRelations->at(i).size() > 0){
				std::vector< Vector3 >  points;
				Coord pa = x[i];
				Coord pb = x[i];
				pb[0] += auxf[i][0];
				pb[1] += auxf[i][1];
				pb[2] += auxf[i][2];
				points.push_back(pa);
				points.push_back(pb);
				Vec<4,float> color2 = Vec<4,float>(1.0f , 1.0f , 0.0f,1.0f);
				vparams->drawTool()->drawLines(points,1,color2 );
				count++;
				Coord e = pa*0.5 + pb*0.5;
				Coord f;
				f[0] = pb[2] - pa[2];
				f[1] = 0;
				f[2] = pa[0] - pb[0];
				std::vector<Vector3> points2;
				points2.push_back(e + f*0.4);
				points2.push_back(pb);
				points2.push_back(e - f*0.4);
				Vec<4,float> color3 = Vec<4,float>(1.0f , 1.0f , 0.0f,1.0f);
				vparams->drawTool()->drawTriangles(points2, color3);
			}
		}
	}
	else{
	}
	

#endif /* SOFA_NO_OPENGL */
}

template <class DataTypes> 
void StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::init()
{
	std::cout << "init FEM" << std::endl;
	this->Inherited::init();
	this->_topology = this->getContext()->getMeshTopology();

        this->tetrahedronInfo.createTopologicalEngine(this->_topology,this->tetrahedronHandler);
        this->tetrahedronInfo.registerTopologicalData();

        this->edgeInfo.createTopologicalEngine(this->_topology);
        this->edgeInfo.registerTopologicalData();

	/** parse the parameter set */
	SetParameterArray paramSet=this->f_parameterSet.getValue();
	if (paramSet.size()>0) {
		this->globalParameters.parameterArray.resize(paramSet.size());
		copy(paramSet.begin(), paramSet.end(),this->globalParameters.parameterArray.begin());
	}
	SetParameterArray paramSet2=f_secondParameterSet.getValue();
	if (paramSet2.size()>0) {
		this->globalParameters2.parameterArray.resize(paramSet2.size());
		copy(paramSet2.begin(), paramSet2.end(),this->globalParameters2.parameterArray.begin());
	}
	/** parse the anisotropy Direction set */
	SetAnisotropyDirectionArray anisotropySet=this->f_anisotropySet.getValue();
	if (anisotropySet.size()>0) {
		this->globalParameters.anisotropyDirection.resize(anisotropySet.size());
		copy(anisotropySet.begin(), anisotropySet.end(),this->globalParameters.anisotropyDirection.begin());
	}
	//vector<HyperelasticMaterialTerm *> materialTermArray;
	/** parse the input material name */
    std::string material = this->f_materialName.getValue();
	if(material=="ArrudaBoyce") {
		fem::BoyceAndArruda<DataTypes> *BoyceAndArrudaMaterial = new fem::BoyceAndArruda<DataTypes>;
		this->myMaterial = BoyceAndArrudaMaterial;
		std::cout<<"The model is "<<material<<std::endl;
		//	materialTermArray =  this->myMaterial->getMaterialTermArray();
	} 
	else if (material=="StVenantKirchhoff"){
		fem::STVenantKirchhoff<DataTypes> *STVenantKirchhoffMaterial = new fem::STVenantKirchhoff<DataTypes>;
		this->myMaterial = STVenantKirchhoffMaterial;
		std::cout<<"The model is "<<material<<std::endl;
		//	materialTermArray =  this->myMaterial->getMaterialTermArray();
	}
	else if (material=="NeoHookean"){
		fem::NeoHookean<DataTypes> *NeoHookeanMaterial = new fem::NeoHookean<DataTypes>;
		this->myMaterial = NeoHookeanMaterial;
		std::cout<<"The model is "<<material<<std::endl;
		//materialTermArray =  this->myMaterial->getMaterialTermArray();
	}
	else if (material=="MooneyRivlin"){
		fem::MooneyRivlin<DataTypes> *MooneyRivlinMaterial = new fem::MooneyRivlin<DataTypes>;
		this->myMaterial = MooneyRivlinMaterial;
		std::cout<<"The model is "<<material<<std::endl;
		//	materialTermArray =  this->myMaterial->getMaterialTermArray();
	}
	else if (material=="VerondaWestman"){
		fem::VerondaWestman<DataTypes> *VerondaWestmanMaterial = new fem::VerondaWestman<DataTypes>;
		this->myMaterial = VerondaWestmanMaterial;
		std::cout<<"The model is "<<material<<std::endl;
		//	materialTermArray =  this->myMaterial->getMaterialTermArray();
	}

	else if (material=="Costa"){
		fem::Costa<DataTypes> *CostaMaterial = new fem::Costa<DataTypes>();
		this->myMaterial =CostaMaterial;
		std::cout<<"The model is "<<material<<std::endl;
		//	materialTermArray =  this->myMaterial->getMaterialTermArray();
	}
	else if (material=="Ogden"){
		fem::Ogden<DataTypes> *OgdenMaterial = new fem::Ogden<DataTypes>();
		this->myMaterial =OgdenMaterial;
		std::cout<<"The model is "<<material<<std::endl;
		//	materialTermArray =  this->myMaterial->getMaterialTermArray();
	}


	else {
        std::cerr << "material name " << material << " is not valid" << std::endl;
	}


	if (!this->_topology->getNbTetrahedra())
	{
        std::cerr << "ERROR(StandardTetrahedralFEMForceField): object must have a Tetrahedral Set Topology." << std::endl;
		return;
	}

	tetrahedronRestInfoVector& tetrahedronInf = *(this->tetrahedronInfo.beginEdit());

	/// prepare to store info in the triangle array
	tetrahedronInf.resize(this->_topology->getNbTetrahedra());

	//helper::vector<typename StandardTetrahedralFEMForceField<DataTypes>::EdgeInformation>& edgeInf = *(this->edgeInfo.beginEdit());
	edgeInformationVector& edgeInf = *(this->edgeInfo.beginEdit());

	edgeInf.resize(this->_topology->getNbEdges());
	this->edgeInfo.endEdit();

	// get restPosition
	if (this->_initialPoints.size() == 0)
	{
		//const VecCoord& p = *this->mstate->getX0();// deprecated ... use: read rest positions
		const VecCoord& p = this->mstate->readRestPositions().ref();
		this->_initialPoints=p;
	}
	int i;

    /// initialize the data structure associate with each tetrahedron
    for (int i=0; i<this->_topology->getNbEdges(); i++)
    {
        edgeInf[i].vertices[0] = (float) this->_topology->getEdge(i)[0];
        edgeInf[i].vertices[1] = (float) this->_topology->getEdge(i)[1];
    }

	/// initialize the data structure associated with each tetrahedron
	for (i=0;i<this->_topology->getNbTetrahedra();++i) {
            this->tetrahedronHandler->applyCreateFunction(i, tetrahedronInf[i],
                        this->_topology->getTetrahedron(i),  (const vector< unsigned int > )0,
                        (const vector< double >)0);
	}
	/// set the call back function upon creation of a tetrahedron

	this->tetrahedronInfo.endEdit();

    /// FOR CUDA
    /// Save the neighbourhood for points (in case of CudaTypes)
    this->initNeighbourhoodPoints();
    this->initNeighbourhoodEdges();

	//testDerivatives();




	
	
	firstime = true; // First loop
	conservationLinearMomentum = true;
	conservationAngularMomentum = true;
	timestep = 0.02;
	counter = -0.02;
	row = -1;

	/*
	// Load internal forces
	 int loadIntForce = -1;
	 internalForce.clear();

	 bool gusano = true;
	 string filename;

	 for(int i=0; i<10; i++)
	 {
		string iSTR = intToString(i);
		if(gusano)
			filename = "C:/Users/amujika/Desktop/CarpetaChachiguay/Gusano simple/Viga/intForcesViga" + iSTR + ".txt";
		else
			filename = "D:/Kodigo/Sofa/share/mesh/forces/IFM1_" + iSTR + ".txt";

		loadIntForce = loadVector(filename, i);

		if(loadIntForce == 0)
			cout << "\nInternal forces " << iSTR << " loaded\n" << endl;
	 }

	
	 */

	m_writeTraces = f_writeTraces.getValue();

	//initialize variables for drawing
	m_gnode = dynamic_cast<simulation::Node*>(this->getContext());

	if(m_gnode){
		sofa::helper::vector<sofa::component::forcefield::LastFriction<defaulttype::Vec3dTypes>*> aux2;
		m_gnode->getTreeObjects<sofa::component::forcefield::LastFriction<defaulttype::Vec3dTypes>>(&aux2);

		if(aux2.size() > 0){
			m_nodes = aux2.at(0)->getNodes();
			m_frictionRelations = aux2.at(0)->getFrictionRelations();
		}
	}

	
	 
} // init()



template <class DataTypes> 
Eigen::Vector3d StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::centerMass(const VecCoord& x, int numVert)
{
	/* Centro de masas calclado tomando todas las masas iguales */
	/*		sum(m_i * r_i)		sum(r_i)
	cm = 	---------------  =	--------
				sum(m_i)		numVert
	
	*/

	Eigen::Matrix<double, 3, 1> cm;
	cm(0) = cm(1) = cm(2) = 0.0;

	for(int i=0; i < x.size(); i++)
	{
		cm(0) += x[i][0];
		cm(1) += x[i][1];
		cm(2) += x[i][2];
	}

	cm(0) = cm(0)/numVert;
	cm(1) = cm(1)/numVert;
	cm(2) = cm(2)/numVert;

	return cm;
}

template <class DataTypes> 
void StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::multLagrange(Eigen::Vector3d &lambdaForce, Eigen::Vector3d &lambdaTorque, const VecCoord& x, Eigen::Vector3d force, Eigen::Vector3d torque, Eigen::Vector3d centerOfMass, int numVert)
{
	Eigen::Matrix3d X, X2, mat, matInv;
	Eigen::Vector3d pos;
	Eigen::Vector3d r, temp, temp2, temp3;
	temp(0) = temp(1) = temp(2) = temp2(0) = temp2(1) = temp2(2) = temp3(0) = temp3(1) = temp3(2) = 0.0;

	for(int i=0; i < x.size(); i++)
	{
		pos(0) = x[i][0];	pos(1) = x[i][1];	pos(2) = x[i][2];
		
		r = pos-centerOfMass;

		temp(0) += r(0); temp(1) += r(1); temp(2) += r(2);
		temp2(0) += r(0)*r(0); temp2(1) += r(1)*r(1); temp2(2) += r(2)*r(2);
		temp3(0) += r(0)*r(1); temp3(1) += r(1)*r(2); temp3(2) += r(0)*r(2);
	}

	X(0,0) = X(1,1) = X(2,2) = 0.0;
	X(0,1) = -temp(2);	X(0,2) = temp(1);
	X(1,0) = temp(2);	X(1,2) = -temp(0);
	X(2,0) = -temp(1);	X(2,1) = temp(0);

	// -X2
	X2(0,0) = temp2(1)+temp2(2);	X2(0,1) = -temp3(0);			X2(0,2) = -temp3(2);
	X2(1,0) = X2(0,1);				X2(1,1) = temp2(0)+temp2(2);	X2(1,2) = -temp3(1);
	X2(2,0) = X2(0,2);				X2(2,1) = X2(1,2);				X2(2,2) = temp2(0)+temp2(1);

	mat = (1.0/numVert) * X * X + X2;
	matInv = mat.inverse();
	
	lambdaTorque = matInv * (-(1.0/numVert) * X * force + torque);

	lambdaForce = (1.0/numVert) * (force + X*lambdaTorque);
}

template <class DataTypes> 
Eigen::Vector3d StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::createForceCorrective(Eigen::Vector3d position, Eigen::Vector3d centerOfMass, Eigen::Vector3d lambdaForce,  Eigen::Vector3d lambdaTorque)
{
	Eigen::Vector3d forceC = -lambdaForce + (position-centerOfMass).cross(lambdaTorque);

	return forceC;
}

template <class DataTypes> 
void StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::addForce(const core::MechanicalParams*  mparams  /* PARAMS FIRST */, DataVecDeriv& d_f, const DataVecCoord& d_x, const DataVecDeriv& /* d_v */)
{

	std::cout << "step: " << getContext()->getTime() << "\n";
    sofa::helper::AdvancedTimer::stepBegin("addForceStandardTetraFEM");

    VecDeriv& f = *d_f.beginEdit();

	auxf = *d_f.beginEdit();

	//std::cout << "before: " << f[52] << "\n";
	//std::cout << "f: " << f[0] << " " << f[1] << " " << f[2] << std::endl;
	
	const VecCoord& x = d_x.getValue();
	
	unsigned int i=0,j=0,k=0,l=0;
	unsigned int nbTetrahedra=this->_topology->getNbTetrahedra();
	
	tetrahedronRestInfoVector& tetrahedronInf = *(this->tetrahedronInfo.beginEdit());


	TetrahedronRestInformation *tetInfo;

	assert(this->mstate);
	this->myposition=x;

	Coord dp[3],x0,sv;
	
	for(i=0; i<nbTetrahedra; i++ )
	{
		tetInfo=&tetrahedronInf[i];
		const Tetrahedron &ta= this->_topology->getTetrahedron(i);
		
		x0=x[ta[0]];
		// compute the deformation gradient
		// deformation gradient = sum of tensor product between vertex position and shape vector
		// optimize by using displacement with first vertex
		dp[0]=x[ta[1]]-x0;
		sv=tetInfo->shapeVector[1];
		for (k=0;k<3;++k) {
			for (l=0;l<3;++l) {
				tetInfo->deformationGradient[k][l]=dp[0][k]*sv[l];
			}
		}
		for (j=1;j<3;++j) {
			dp[j]=x[ta[j+1]]-x0;
			sv=tetInfo->shapeVector[j+1];
			for (k=0;k<3;++k) {
				for (l=0;l<3;++l) {
					tetInfo->deformationGradient[k][l]+=dp[j][k]*sv[l];
				}
			}
		}
		
		//Compute the matrix strain displacement B 6*3
		for (int alpha=0; alpha<4; ++alpha){
			Coord sva=tetInfo->shapeVector[alpha];
			Matrix63 matBa;
			matBa[0][0]=tetInfo->deformationGradient[0][0]*sva[0];
			matBa[0][1]=tetInfo->deformationGradient[1][0]*sva[0];
			matBa[0][2]=tetInfo->deformationGradient[2][0]*sva[0];

			matBa[2][0]=tetInfo->deformationGradient[0][1]*sva[1];
			matBa[2][1]=tetInfo->deformationGradient[1][1]*sva[1];
			matBa[2][2]=tetInfo->deformationGradient[2][1]*sva[1];

			matBa[5][0]=tetInfo->deformationGradient[0][2]*sva[2];
			matBa[5][1]=tetInfo->deformationGradient[1][2]*sva[2];
			matBa[5][2]=tetInfo->deformationGradient[2][2]*sva[2];

			matBa[1][0]=(tetInfo->deformationGradient[0][0]*sva[1]+tetInfo->deformationGradient[0][1]*sva[0]);
			matBa[1][1]=(tetInfo->deformationGradient[1][0]*sva[1]+tetInfo->deformationGradient[1][1]*sva[0]);
			matBa[1][2]=(tetInfo->deformationGradient[2][0]*sva[1]+tetInfo->deformationGradient[2][1]*sva[0]);

			matBa[3][0]=(tetInfo->deformationGradient[0][2]*sva[0]+tetInfo->deformationGradient[0][0]*sva[2]);
			matBa[3][1]=(tetInfo->deformationGradient[1][2]*sva[0]+tetInfo->deformationGradient[1][0]*sva[2]);
			matBa[3][2]=(tetInfo->deformationGradient[2][2]*sva[0]+tetInfo->deformationGradient[2][0]*sva[2]);

			matBa[4][0]=(tetInfo->deformationGradient[0][1]*sva[2]+tetInfo->deformationGradient[0][2]*sva[1]);
			matBa[4][1]=(tetInfo->deformationGradient[1][1]*sva[2]+tetInfo->deformationGradient[1][2]*sva[1]);
			matBa[4][2]=(tetInfo->deformationGradient[2][1]*sva[2]+tetInfo->deformationGradient[2][2]*sva[1]);

			tetInfo->matB[alpha]=matBa;
		}


		/// compute the right Cauchy-Green deformation matrix
		for (k=0;k<3;++k) {
			for (l=k;l<3;++l) {
				//	if (l>=k) {
				tetInfo->deformationTensor(k,l)=(tetInfo->deformationGradient(0,k)*tetInfo->deformationGradient(0,l)+
					tetInfo->deformationGradient(1,k)*tetInfo->deformationGradient(1,l)+
					tetInfo->deformationGradient(2,k)*tetInfo->deformationGradient(2,l));
				//	}
				//	else 
				//		tetInfo->deformationTensor[k][l]=tetInfo->deformationTensor[l][k];
			}
		}


		//in case of transversaly isotropy

		if(this->globalParameters.anisotropyDirection.size()>0){
			tetInfo->fiberDirection=this->globalParameters.anisotropyDirection[0];
			Coord vectCa=tetInfo->deformationTensor*tetInfo->fiberDirection;
			Real aDotCDota=dot(tetInfo->fiberDirection,vectCa);
			tetInfo->lambda=(Real)sqrt(aDotCDota);
		}
		Coord areaVec = cross( dp[1], dp[2] );

		tetInfo->J = dot( areaVec, dp[0] ) * tetInfo->volScale;
		tetInfo->trC = (Real)( tetInfo->deformationTensor(0,0) + tetInfo->deformationTensor(1,1) + tetInfo->deformationTensor(2,2));
		//tetInfo->strainEnergy=this->myMaterial->getStrainEnergy(tetInfo,this->globalParameters); // to uncomment for test derivatives
		tetInfo->SPKTensorGeneral.clear();
		MatrixSym SPK;
		if(i < f_secondMaterialStart.getValue()){
			this->myMaterial->deriveSPKTensor(tetInfo,this->globalParameters,SPK); // calculate the SPK tensor of the chosen material
		}
		else{
			this->myMaterial->deriveSPKTensor(tetInfo,this->globalParameters2,SPK); 
		}
		tetInfo->SPKTensorGeneral=SPK;
		if(i==0){
			//cout << "prueba: " << sv << "\n" << tetInfo->matB[0].transposed()  <<  "\n" << SPK << "\n" << tetInfo->restVolume << std::endl;
		}
		for(l=0;l<4;++l){
			f[ta[l]]-=tetInfo->matB[l].transposed()*SPK*tetInfo->restVolume;
		}
	}

	for(unsigned int i = 0; i < f.size(); i++){
		auxf[i] += f[i];
		if(i == 10){
			std::cout << "material: " << auxf[i] << "\n";
		}
		float module = sqrt(auxf[i][0]*auxf[i][0] + auxf[i][1]*auxf[i][1] + auxf[i][2]*auxf[i][2]);
		auxf[i][0] = auxf[i][0]/module;
		auxf[i][1] = auxf[i][1]/module;
		auxf[i][2] = auxf[i][2]/module;
	}
	

	/// indicates that the next call to addDForce will need to update the stiffness matrix
	this->updateMatrix=true;
	this->tetrahedronInfo.endEdit();
	//std::cout << "f0: " << f[0] << " " << f[1] << " " << f[2] << std::endl;

	d_f.endEdit();

    sofa::helper::AdvancedTimer::stepEnd("addForceStandardTetraFEM");

	//std::cout << "standard: " << f[52] << "\n";

	
}


template <class DataTypes> 
void StandardTetrahedralFEMForceFieldTwoMaterials<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& d_df, const DataVecDeriv& d_dx)
{
	
	if(m_writeTraces) std::cout << "standard" << std::endl;
	/*double incremento = 0.02;

	antTimestep = timestep;

	if(firstime)
	{
		timestep = 0.0;
		firstime = false;
	}
	else
		timestep += incremento;


	if(timestep+incremento == PERIOD || PERIOD < timestep+incremento)
		firstime = true; // reseteamos el temistep para que vuelva a tomar la primera fila

	//cout << "\ntimestep: " << timestep << endl;
	
	row = floor(8.0f * timestep / PERIOD);

	vector<double> interpolateForces = interpolate(row);*/
	

    sofa::helper::AdvancedTimer::stepBegin("addDForceStandardTetraFEM");

    VecDeriv& df = *d_df.beginEdit();

	//cout << "\nSara df addDD " << df[0][0] << " " << df[0][1] << " " << df[0][2] << endl;

	const VecDeriv& dx = d_dx.getValue();
	Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

	//cout << "\nSara df addDDinter " << df[0][0] << " " << df[0][1] << " " << df[0][2] << endl;

	unsigned int i=0,j=0,k=0,l=0;
	unsigned int nbEdges=this->_topology->getNbEdges();
	
	const vector< Edge> &edgeArray=this->_topology->getEdges() ;

	helper::vector<EdgeInformation>& edgeInf = *(this->edgeInfo.beginEdit());
	tetrahedronRestInfoVector& tetrahedronInf = *(this->tetrahedronInfo.beginEdit());

	EdgeInformation *einfo;


	/// if the  matrix needs to be updated
	if (this->updateMatrix) {

		TetrahedronRestInformation *tetInfo;
		unsigned int nbTetrahedra=this->_topology->getNbTetrahedra();
		const std::vector< Tetrahedron> &tetrahedronArray=this->_topology->getTetrahedra() ;


		for(l=0; l<nbEdges; l++ )edgeInf[l].DfDx.clear();
		for(i=0; i<nbTetrahedra; i++ )
		{
			tetInfo=&tetrahedronInf[i];
//			Matrix3 &df=tetInfo->deformationGradient;
//			Matrix3 Tdf=df.transposed();
			BaseMeshTopology::EdgesInTetrahedron te=this->_topology->getEdgesInTetrahedron(i);

			/// describe the jth vertex index of triangle no i 
			const Tetrahedron &ta= tetrahedronArray[i];
			for(j=0;j<6;j++) {
				einfo= &edgeInf[te[j]];
				Edge e=this->_topology->getLocalEdgesInTetrahedron(j);

				k=e[0];
				l=e[1];
				if (edgeArray[te[j]][0]!=ta[k]) {
					k=e[1];
					l=e[0];
				}
				Matrix3 &edgeDfDx = einfo->DfDx;
				
				Coord svl=tetInfo->shapeVector[l];
				Coord svk=tetInfo->shapeVector[k];

				Matrix3  M, N;
				Matrix6 outputTensor;
				N.clear();

				// Calculates the dS/dC tensor 6*6
				if(i < f_secondMaterialStart.getValue()){
					this->myMaterial->ElasticityTensor(tetInfo,this->globalParameters,outputTensor); // calculate the SPK tensor of the chosen material
				}
				else{
					this->myMaterial->ElasticityTensor(tetInfo,this->globalParameters2,outputTensor);
				}
				Matrix63 mBl=tetInfo->matB[l];
				mBl[1][0]/=2;mBl[1][1]/=2;mBl[1][2]/=2;mBl[3][0]/=2;mBl[3][1]/=2;mBl[3][2]/=2;mBl[4][0]/=2;mBl[4][1]/=2;mBl[4][2]/=2;

				N=(tetInfo->matB[k].transposed()*outputTensor*mBl);




				//Now M
				Real productSD=0;

				Coord vectSD=tetInfo->SPKTensorGeneral*svk;
				productSD=dot(vectSD,svl);
				M[0][1]=M[0][2]=M[1][0]=M[1][2]=M[2][0]=M[2][1]=0;
				M[0][0]=M[1][1]=M[2][2]=(Real)productSD;

				edgeDfDx += (M+N.transposed())*tetInfo->restVolume;


			}// end of for j				
		}//end of for i
		this->updateMatrix=false;
	}// end of if


	/// performs matrix vector computation
	unsigned int v0,v1;
	Deriv deltax;	Deriv dv0,dv1;

	for(l=0; l<nbEdges; l++ )
	{
		einfo=&edgeInf[l];
		v0=edgeArray[l][0];
		v1=edgeArray[l][1];

		deltax= dx[v0] - dx[v1];
		dv0 = einfo->DfDx * deltax;
		// do the transpose multiply:
		dv1[0] = (Real)(deltax[0]*einfo->DfDx[0][0] + deltax[1]*einfo->DfDx[1][0] + deltax[2]*einfo->DfDx[2][0]);
		dv1[1] = (Real)(deltax[0]*einfo->DfDx[0][1] + deltax[1]*einfo->DfDx[1][1] + deltax[2]*einfo->DfDx[2][1]);
		dv1[2] = (Real)(deltax[0]*einfo->DfDx[0][2] + deltax[1]*einfo->DfDx[1][2] + deltax[2]*einfo->DfDx[2][2]);
		// add forces
		df[v0] += dv1 * kFactor;
		df[v1] -= dv0 * kFactor;

	}

	//cout << "\nSara df final addDD " << df[0][0] << " " << df[0][1] << " " << df[0][2] << endl;

	this->edgeInfo.endEdit();
	this->tetrahedronInfo.endEdit();
	d_df.beginEdit();

    sofa::helper::AdvancedTimer::stepEnd("addDForceStandardTetraFEM");
}


} // namespace forcefield
} // namespace component
} // namespace sofa

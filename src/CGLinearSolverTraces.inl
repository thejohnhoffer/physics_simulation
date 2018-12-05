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
#ifndef SOFA_COMPONENT_LINEARSOLVER_CGLinearSolverTraces_INL
#define SOFA_COMPONENT_LINEARSOLVER_CGLinearSolverTraces_INL

#include "CGLinearSolverTraces.h"
#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/SparseMatrix.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.h>
#include <sofa/simulation/common/MechanicalVisitor.h>
#include <sofa/helper/system/thread/CTime.h>
#include <sofa/helper/AdvancedTimer.h>

#include <sofa/core/ObjectFactory.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace linearsolver
{

/// Linear system solver using the conjugate gradient iterative algorithm
template<class TMatrix, class TVector>
CGLinearSolverTraces<TMatrix,TVector>::CGLinearSolverTraces()
	:  f_writeTraces( initData(&f_writeTraces, true, "write_traces","write traces in console") )
{

}

/// Solve Mx=b
template<class TMatrix, class TVector>
void CGLinearSolverTraces<TMatrix,TVector>::solve(Matrix& M, Vector& x, Vector& b)
{
	//sofa::component::linearsolver::GraphScatteredMatrix myM = (sofa::component::linearsolver::GraphScatteredMatrix) (M);
	std::cout << "aber: " << typeid(M).name() << "\n";
	std::cout << "aber2: " << &M << "\n";
	bool writeTraces = f_writeTraces.getValue();
#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printComment("ConjugateGradient");
#endif

#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printNode("VectorAllocation");
#endif
    const core::ExecParams* params = core::ExecParams::defaultInstance();
    typename Inherit::TempVectorContainer vtmp(this, params, M, x, b);
    Vector& p = *vtmp.createTempVector();
    Vector& q = *vtmp.createTempVector();
    Vector& r = *vtmp.createTempVector();

    const bool printLog = this->f_printLog.getValue();
    const bool verbose  = Inherit::f_verbose.getValue();

    // -- solve the system using a conjugate gradient solution
    double rho, rho_1=0, alpha, beta;

    if( verbose )
        sout<<"CGLinearSolver, b = "<< b <<sendl;

    if( Inherit::f_warmStart.getValue() )
    {
        r = M * x;
        r.eq( b, r, -1.0 );   //  initial residual r = b - Ax;
    }
    else
    {
        x.clear();
        r = b; // initial residual
    }

    double normb = b.norm();
    std::map < std::string, sofa::helper::vector<SReal> >& graph = *Inherit::f_graph.beginEdit();
    sofa::helper::vector<SReal>& graph_error = graph[(this->isMultiGroup()) ? this->currentNode->getName()+std::string("-Error") : std::string("Error")];
    graph_error.clear();
    sofa::helper::vector<SReal>& graph_den = graph[(this->isMultiGroup()) ? this->currentNode->getName()+std::string("-Denominator") : std::string("Denominator")];
    graph_den.clear();
    graph_error.push_back(1);
    unsigned nb_iter;
    const char* endcond = "iterations";

#ifdef DISPLAY_TIME
    sofa::helper::system::thread::CTime timer;
    time1 = (SReal) timer.getTime();
#endif

#ifdef SOFA_DUMP_VISITOR_INFO
    simulation::Visitor::printCloseNode("VectorAllocation");
#endif
    for( nb_iter=1; nb_iter<=Inherit::f_maxIter.getValue(); nb_iter++ )
    {
		//if(writeTraces) std::cout << "iter: " << nb_iter << std::endl;
      //std::cout << "Iteration_" << nb_iter;
#ifdef SOFA_DUMP_VISITOR_INFO
        std::ostringstream comment;
        if (simulation::Visitor::isPrintActivated())
        {
            comment << "Iteration_" << nb_iter;
            simulation::Visitor::printNode(comment.str());
        }
#endif
        // 		printWithElapsedTime( x, helper::system::thread::CTime::getTime()-time0,sout );

        //z = r; // no precond
        //rho = r.dot(z);
        rho = r.dot(r);

        if (nb_iter>1)
        {
            double normr = sqrt(rho); //sqrt(r.dot(r));
            double err = normr/normb;
            graph_error.push_back(err);
			if(writeTraces) std::cout << "tol: " << err << " " << Inherit::f_tolerance.getValue() << std::endl;
            if (err <= Inherit::f_tolerance.getValue())
            {
                endcond = "tolerance";

#ifdef SOFA_DUMP_VISITOR_INFO
                if (simulation::Visitor::isPrintActivated())
                    simulation::Visitor::printCloseNode(comment.str());
#endif
                break;
            }
        }

        if( nb_iter==1 )
            p = r; //z;
        else
        {
            beta = rho / rho_1;
            //p = p*beta + r; //z;
            //Inherit::cgstep_beta(params, p,r,beta);
	    p *= beta;
	    p += r; //z;
        }

        if( verbose )
        {
            sout<<"p : "<<p<<sendl;
        }
        // matrix-vector product
        q = M*p;

        if( verbose )
        {
            sout<<"q = M p : "<<q<<sendl;
        }

        double den = p.dot(q);

        graph_den.push_back(den);

		//if(writeTraces) std::cout << "den: " << fabs(den) << " " << Inherit::f_smallDenominatorThreshold.getValue() << std::endl;

        if( fabs(den)<Inherit::f_smallDenominatorThreshold.getValue() )
        {
            endcond = "threshold";
            if( verbose )
            {
                sout<<"CGLinearSolver, den = "<<den<<", smallDenominatorThreshold = "<<Inherit::f_smallDenominatorThreshold.getValue()<<sendl;
            }
#ifdef SOFA_DUMP_VISITOR_INFO
            if (simulation::Visitor::isPrintActivated())
                simulation::Visitor::printCloseNode(comment.str());
#endif
            break;
        }
        alpha = rho/den;
        //x.peq(p,alpha);                 // x = x + alpha p
        //r.peq(q,-alpha);                // r = r - alpha q
        //Inherit::cgstep_alpha(params, x,r,p,q,alpha);
	x.peq(p,alpha);                 // x = x + alpha p
	r.peq(q,-alpha);                // r = r - alpha q
        if( verbose )
        {
            sout<<"den = "<<den<<", alpha = "<<alpha<<sendl;
            sout<<"x : "<<x<<sendl;
            sout<<"r : "<<r<<sendl;
        }

        rho_1 = rho;
#ifdef SOFA_DUMP_VISITOR_INFO
        if (simulation::Visitor::isPrintActivated())
            simulation::Visitor::printCloseNode(comment.str());
#endif
    }

#ifdef DISPLAY_TIME
    time1 = (SReal)(((SReal) timer.getTime() - time1) * timeStamp / (nb_iter-1));
#endif

    Inherit::f_graph.endEdit();

    sofa::helper::AdvancedTimer::valSet("CG iterations", nb_iter);

    // x is the solution of the system
    if( printLog )
    {
#ifdef DISPLAY_TIME
        std::cerr<<"CGLinearSolver::solve, CG = "<<time1<<" build = "<<time2<<std::endl;
#endif
        sout<<"CGLinearSolver::solve, nbiter = "<<nb_iter<<" stop because of "<<endcond<<sendl;
    }
    if( verbose )
    {
        sout<<"CGLinearSolver::solve, solution = "<<x<<sendl;
    }
    vtmp.deleteTempVector(&p);
    vtmp.deleteTempVector(&q);
    vtmp.deleteTempVector(&r);
}


} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_LINEARSOLVER_CGLINEARSOLVER_INL

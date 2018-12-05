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
#ifndef SOFA_COMPONENT_LINEARSOLVER_CGLinearSolverTraces_H
#define SOFA_COMPONENT_LINEARSOLVER_CGLinearSolverTraces_H

#include <SofaBaseLinearSolver/CGLinearSolver.h>

#include <sofa/helper/map.h>

#include <math.h>

#include "initPlugin.h"


namespace sofa
{

namespace component
{

namespace linearsolver
{

//#define DISPLAY_TIME

/// Linear system solver using the conjugate gradient iterative algorithm
template<class TMatrix, class TVector>
class CGLinearSolverTraces : public sofa::component::linearsolver::CGLinearSolver<TMatrix, TVector>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(CGLinearSolverTraces,TMatrix,TVector),SOFA_TEMPLATE2(sofa::component::linearsolver::CGLinearSolver,TMatrix,TVector));

    typedef CGLinearSolver<TMatrix, TVector> Inherit;
    typedef TMatrix Matrix;
    typedef TVector Vector;
#ifdef DISPLAY_TIME
    SReal time1;
    SReal time2;
    SReal timeStamp;
#endif
protected:

    CGLinearSolverTraces();

public:

    /// Solve Mx=b
    void solve (Matrix& M, Vector& x, Vector& b);

	Data<bool> f_writeTraces;

};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_BUILD_SIELEGANSPLUGIN)
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< GraphScatteredMatrix, GraphScatteredVector >;
#ifndef SOFA_FLOAT
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< FullMatrix<double>, FullVector<double> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< SparseMatrix<double>, FullVector<double> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<double>, FullVector<double> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<2,2,double> >, FullVector<double> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<3,3,double> >, FullVector<double> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<4,4,double> >, FullVector<double> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<6,6,double> >, FullVector<double> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<8,8,double> >, FullVector<double> >;
#endif

#ifndef SOFA_DOUBLE
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<float>, FullVector<float> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<2,2,float> >, FullVector<float> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<3,3,float> >, FullVector<float> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<4,4,float> >, FullVector<float> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<6,6,float> >, FullVector<float> >;
extern template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<8,8,float> >, FullVector<float> >;
#endif
#endif

} // namespace linearsolver

} // namespace component

} // namespace sofa

#endif

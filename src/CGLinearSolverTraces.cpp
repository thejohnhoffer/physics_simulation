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
// Author: François Faure, INRIA-UJF, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
#include "CGLinearSolverTraces.inl"

#include <sofa/core/ObjectFactory.h>
#include <iostream>
#include <SofaBaseLinearSolver/FullMatrix.h>
#include <SofaBaseLinearSolver/FullVector.h>
#include <SofaBaseLinearSolver/GraphScatteredTypes.h>

namespace sofa
{

namespace component
{

namespace linearsolver
{

using namespace sofa::defaulttype;
using sofa::core::MultiVecDerivId;


SOFA_DECL_CLASS(CGLinearSolverTraces)

int CGLinearSolverTracesClass = core::RegisterObject("Linear system solver using the conjugate gradient iterative algorithm")
        .add< CGLinearSolverTraces< GraphScatteredMatrix, GraphScatteredVector > >(true)
#ifndef SOFA_FLOAT
        .add< CGLinearSolverTraces< FullMatrix<double>, FullVector<double> > >()
        .add< CGLinearSolverTraces< SparseMatrix<double>, FullVector<double> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<double>, FullVector<double> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<2,2,double> >, FullVector<double> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<3,3,double> >, FullVector<double> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<4,4,double> >, FullVector<double> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<6,6,double> >, FullVector<double> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<8,8,double> >, FullVector<double> > >()
#endif
#ifndef SOFA_DOUBLE
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<float>, FullVector<float> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<2,2,float> >, FullVector<float> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<3,3,float> >, FullVector<float> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<4,4,float> >, FullVector<float> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<6,6,float> >, FullVector<float> > >()
        .add< CGLinearSolverTraces< CompressedRowSparseMatrix<Mat<8,8,float> >, FullVector<float> > >()
#endif
        .addAlias("MyCGSolver")
        .addAlias("MyConjugateGradient")
        ;

template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< GraphScatteredMatrix, GraphScatteredVector >;
#ifndef SOFA_FLOAT
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< FullMatrix<double>, FullVector<double> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< SparseMatrix<double>, FullVector<double> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<double>, FullVector<double> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<2,2,double> >, FullVector<double> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<3,3,double> >, FullVector<double> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<4,4,double> >, FullVector<double> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<6,6,double> >, FullVector<double> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<8,8,double> >, FullVector<double> >;
#endif

#ifndef SOFA_DOUBLE
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<float>, FullVector<float> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<2,2,float> >, FullVector<float> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<3,3,float> >, FullVector<float> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<4,4,float> >, FullVector<float> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<6,6,float> >, FullVector<float> >;
template class SOFA_SiElegansPlugin_API CGLinearSolverTraces< CompressedRowSparseMatrix<defaulttype::Mat<8,8,float> >, FullVector<float> >;
#endif
} // namespace linearsolver

} // namespace component

} // namespace sofa


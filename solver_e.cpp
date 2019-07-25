/*
* Copyright (c) 2017, 
* Eyal Shalev (eyal@gsi.gov.il)
* Vladimir Lyakhovsky
* Harel Levin (harellevin@gmail.com)
* Gal Oren (galoren.com@gmail.com)
* All rights reserved to:
* Geological Survey of Israel (GSI) &
* Nuclear Research Center - Negev (NRCN).
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*    * Redistributions of source code must retain the above copyright
*      notice, this list of conditions and the following disclaimer.
*    * Redistributions in binary form must reproduce the above copyright
*      notice, this list of conditions and the following disclaimer in the
*      documentation and/or other materials provided with the distribution.
*    * Neither the name of Harel Levin or Gal Oren, nor the
*      names of its contributors may be used to endorse or promote products
*      derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Harel Levin & Gal Oren BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#define _GLIBCXX_USE_CXX11_ABI 0

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_SerialComm.h>
#include <Epetra_CrsMatrix.h>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

using std::cout;
using std::endl;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::ParameterList;
using Teuchos::parameterList;

typedef Epetra_Map map_type;
typedef double scalar_type;
typedef int local_ordinal_type;
typedef int global_ordinal_type;
typedef Epetra_CrsMatrix crs_matrix_type;
typedef Epetra_CrsGraph crs_graph_type;
typedef Epetra_Vector vector_type;
typedef Epetra_MultiVector multivector_type;
typedef Epetra_Operator operator_type;

RCP<const Epetra_Comm > comm;
ArrayRCP<size_t> rowPtr;
RCP<crs_matrix_type> a_matrix;
ArrayRCP<int> numOfIndicesPerRow;
ArrayRCP<global_ordinal_type> colInd;
int numGlobalEntries;
int numOfNnz;

extern "C"{

void initialize_solver_ (int *a_n, int *a_nnz, size_t *a_rowptr, int *a_colind);
void run_solver_gmres_ (double *a_values, double *b_values, int *info);
void run_solver_cg_ (double *a_values, double *b_values, int *info);
void run_solver_ (std::string solverType, double *a_values, double *b_values, int *info);
}

void initialize_solver_(int *a_n, int *a_nnz, size_t *a_rowptr, int *a_colind)
{
  comm  = RCP<const Epetra_Comm >(new Epetra_SerialComm);

  global_ordinal_type indexBase = 0;
  numGlobalEntries = *a_n;
  numOfNnz = *a_nnz;
  const RCP<const map_type> rowMap =
    rcp (new map_type (numGlobalEntries, indexBase, *comm));
  const RCP<const map_type> colMap =
    rcp (new map_type (numGlobalEntries, indexBase, *comm));

  rowPtr = ArrayRCP<size_t>(numGlobalEntries+1);
  rowPtr.deepCopy(Teuchos::ArrayView<size_t>(a_rowptr, numGlobalEntries+1));
  colInd =  ArrayRCP<global_ordinal_type>(numOfNnz);
  colInd.deepCopy(Teuchos::ArrayView<global_ordinal_type>(a_colind, numOfNnz));
  
  numOfIndicesPerRow = ArrayRCP<int>(numGlobalEntries);
  //Change Fortran indexing to CPP
  #pragma omp parallel
  {
    int i;
#if !defined(_AIX)
    #pragma omp for simd
#else
    #pragma omp for
#endif
    for (i=0;i<=numGlobalEntries;++i)
    {
      rowPtr[i]--;
    }
#if !defined(_AIX)
    #pragma omp for simd
#else
    #pragma omp for
#endif
    for (i=1;i<=numGlobalEntries;++i)
    {
      numOfIndicesPerRow[i-1] = rowPtr[i] - rowPtr[i-1];
    }
#if !defined(_AIX)
    #pragma omp for simd
#else
    #pragma omp for
#endif
    for (i=0;i<numOfNnz;++i)
    {
      colInd[i]--;
    }
  }
  

  RCP<crs_graph_type> crsGraph = rcp(new crs_graph_type (Copy, *rowMap, numOfIndicesPerRow.getRawPtr(), true));

  const RCP<ParameterList> crsGraphParameters =
    rcp (new ParameterList("CRSGraphParams"));

#if !defined(_AIX)
    #pragma omp parallel for simd
#else
    #pragma omp parallel for
#endif
  for (int i=0;i<numGlobalEntries;++i)
  {
    int status = crsGraph->InsertGlobalIndices(i, numOfIndicesPerRow[i], colInd.getRawPtr() + rowPtr[i]);
  }

  crsGraph->FillComplete();

  a_matrix=RCP<crs_matrix_type>(new crs_matrix_type(Copy, *crsGraph));
}

void run_solver_gmres_(double *a_values, double *b_values, int *info)
{
    run_solver_("GMRES", a_values, b_values, info);
}

void run_solver_cg_(double *a_values, double *b_values, int *info)
{
    run_solver_("CG", a_values, b_values, info);
}

void run_solver_(std::string solverType, double *a_values, double *b_values, int *info)
{
    const global_ordinal_type indexBase = 0;

    RCP<const map_type> rowMap =
      rcp (new map_type (numGlobalEntries, indexBase, *comm));
    RCP<const map_type> colMap =
      rcp (new map_type (numGlobalEntries, indexBase, *comm));

    ArrayRCP<scalar_type> values (a_values, 0, numOfNnz, false);

    for(size_t i=indexBase;i<indexBase+numGlobalEntries;i++)
    {
      a_matrix->ReplaceGlobalValues(i, rowPtr[i - indexBase +1]-rowPtr[i - indexBase],a_values+rowPtr[i - indexBase],colInd.getRawPtr()+rowPtr[i - indexBase]);
    }
    a_matrix->FillComplete();
    
    Teuchos::ArrayView<const scalar_type> b_vector(b_values, numGlobalEntries);
    Teuchos::ArrayView<const Teuchos::ArrayView<const scalar_type>> b_arrayview(&b_vector, 1);
    double *b_array_ptr = const_cast<double*>(b_arrayview[0].getRawPtr());
    RCP<multivector_type> b_multivector = 
        rcp (new multivector_type(Copy, *rowMap, &b_array_ptr, 1));

    RCP<multivector_type> x_multivector = 
        rcp (new multivector_type(Copy, a_matrix->DomainMap(), &b_values, 1));

    RCP<ParameterList> solverParams = parameterList();
    solverParams->set ("Num Blocks", 40);
    //solverParams->set ("Recycled Blocks", 20);
    solverParams->set ("Maximum Iterations", 1000);
    solverParams->set ("Convergence Tolerance", 1.0e-5);

    Belos::SolverFactory<scalar_type, multivector_type, operator_type> factory;
    RCP<Belos::SolverManager<scalar_type, multivector_type, operator_type>> solver = 
        factory.create(solverType, solverParams);
    
    typedef Belos::LinearProblem<scalar_type, multivector_type, operator_type> problem_type;
    RCP<problem_type> problem = rcp (new problem_type(a_matrix, x_multivector, b_multivector));
    problem->setProblem();

    solver->setProblem(problem);

    Belos::ReturnType result = solver->solve();
    #pragma omp barrier
    const int numIters = solver->getNumIters();
    if(result == Belos::Converged)
    {
        std::cout << "CONVERGED! :) in " << numIters << " iterations.\n";
    }
    else
    {
        std::cout << "DIDN'T CONVERGED :(\n";
    }

    ArrayRCP<const scalar_type> new_b ((*x_multivector)[0], 0, numGlobalEntries, false);
#if !defined(_AIX)
    #pragma omp parallel for simd
#else
    #pragma omp parallel for
#endif
    for (int i=0; i < numGlobalEntries; ++i)
    {
        b_values[i] = new_b[i];
    }
}

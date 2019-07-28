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
*    * Neither the name of Eyal Shalev, Vladimir Lyakhovsky, Harel Levin or 
*      Gal Oren, nor the names of its contributors may be used to endorse 
*      or promote products derived from this software without specific prior 
*      written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Eyal Shalev, Vladimir Lyakhovsky, Harel Levin 
* & Gal Oren BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  ON ANY THEORY OF LIABILITY, 
* WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
* OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
* ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>

using std::cout;
using std::endl;
using Teuchos::Comm;
using Teuchos::SerialComm;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::ParameterList;
using Teuchos::parameterList;

typedef Tpetra::Map<> map_type;
typedef Tpetra::Vector<>::scalar_type scalar_type;
typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;
typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;
typedef Tpetra::CrsMatrix<> crs_matrix_type;
typedef Tpetra::CrsGraph<> crs_graph_type;
typedef Tpetra::Vector<> vector_type;
typedef Tpetra::MultiVector<> multivector_type;
typedef Tpetra::Operator<> operator_type;

RCP<const Comm<int> > comm;
RCP<crs_matrix_type> a_matrix;
ArrayRCP<size_t> rowPtr;
ArrayRCP<global_ordinal_type> colInd;
Tpetra::global_size_t numGlobalEntries;
Tpetra::global_size_t numOfNnz;

extern "C"{

void initialize_solver_ (int *a_n, int *a_nnz, size_t *a_rowptr, int *a_colind);
void run_solver_gmres_ (double *a_values, double *b_values, int *info);
void run_solver_cg_ (double *a_values, double *b_values, int *info);
void run_solver_ (std::string solverType, double *a_values, double *b_values, int *info);
void finalize_solver_ ();

}

void initialize_solver_(int *a_n, int *a_nnz, size_t *a_rowptr, int *a_colind)
{
  comm  = RCP<const Comm<int> >(new SerialComm<int> ());

  const global_ordinal_type indexBase = 0;
  numGlobalEntries = *a_n;
  numOfNnz = *a_nnz;
  const RCP<const map_type> rowMap =
    rcp (new map_type (numGlobalEntries, indexBase, comm));
  const RCP<const map_type> colMap =
    rcp (new map_type (numGlobalEntries, indexBase, comm));

  rowPtr = ArrayRCP<size_t>(numGlobalEntries+1);
  rowPtr.deepCopy(Teuchos::ArrayView<size_t>(a_rowptr, numGlobalEntries+1));
  colInd =  ArrayRCP<global_ordinal_type>(numOfNnz);
  colInd.deepCopy(Teuchos::ArrayView<global_ordinal_type>(a_colind, numOfNnz));
  
  //Change Fortran indexing to CPP
  #pragma omp parallel
  {
    int i;
    #pragma omp for simd
    for (i=0;i<numGlobalEntries;++i)
    {
      rowPtr[i]--;
    }
    #pragma omp for simd
    for (i=0;i<numOfNnz;++i)
    {
      colInd[i]--;
    }
  }
  

  RCP<crs_graph_type> crsGraph = 
    rcp (new crs_graph_type (rowMap, colMap, rowPtr, colInd));

  const RCP<ParameterList> crsGraphParameters =
    rcp (new ParameterList("CRSGraphParams"));

  crsGraph->fillComplete();

  a_matrix=RCP<crs_matrix_type>(new crs_matrix_type(crsGraph));
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
      rcp (new map_type (numGlobalEntries, indexBase, comm));
    RCP<const map_type> colMap =
      rcp (new map_type (numGlobalEntries, indexBase, comm));

    ArrayRCP<scalar_type> values (a_values, 0, numOfNnz, false);

    a_matrix->resumeFill();
    for(size_t i=indexBase;i<indexBase+numGlobalEntries;i++)
    {
      a_matrix->replaceGlobalValues(i, rowPtr[i - indexBase +1]-rowPtr[i - indexBase],a_values+rowPtr[i - indexBase],colInd.getRawPtr()+rowPtr[i - indexBase]);
    }
    a_matrix->fillComplete();
    a_matrix->print(std::cout);

    RCP<multivector_type> x_multivector = 
        rcp (new multivector_type(a_matrix->getDomainMap(), 1));
    x_multivector->print(std::cout);
    
    Teuchos::ArrayView<const scalar_type> b_vector(b_values, numGlobalEntries);
    Teuchos::ArrayView<const Teuchos::ArrayView<const scalar_type>> b_arrayview(&b_vector, 1);
    RCP<multivector_type> b_multivector = 
        rcp (new multivector_type(a_matrix->getRangeMap(), b_arrayview, 1));
    b_multivector->print(std::cout);

    deep_copy(*x_multivector, *b_multivector);

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

    ArrayRCP<const scalar_type> new_b = x_multivector->get1dView();
    #pragma omp parallel for simd
    for (int i=0; i < numGlobalEntries; ++i)
    {
        b_values[i] = new_b[i];
    }
}

void finalize_solver_()
{
    a_matrix.release();
//    a_matrix=RCP<crs_matrix_type>(new crs_matrix_type(Teuchos::null));
}

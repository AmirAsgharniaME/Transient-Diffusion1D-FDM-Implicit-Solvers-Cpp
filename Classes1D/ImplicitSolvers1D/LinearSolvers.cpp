#include "ImplicitSolvers1D/LinearSolvers.hpp"

void LinearSolvers::Solve_Gaussian_Elimination(CoefficientMatrix& A, RHS& RHS_Obj, Field1D& Solution)

 {
    const std::size_t NumRows = A.GetNumRows();
    const std::size_t NumCols = A.GetNumColumns();

    if (NumRows != NumCols)
    {
        throw std::invalid_argument("Solve: matrix A must be square.");
    }

    if (RHS_Obj.GetSize() != NumRows || Solution.Size() != NumRows)
    {
        throw std::invalid_argument("Solve: size of RHS and Solution must match matrix size.");
    }

    GaussianElimination::ForwardElimination(A, RHS_Obj);
    GaussianElimination::BackSubstitution(A, RHS_Obj, Solution);

}


void LinearSolvers::Solve_Thomas_Algorithm(TridiagonalMatrix& A3_Obj, RHS& RHS_Obj, Field1D& Solution)
{
   ThomasAlgorithm::ForwardElimination(A3_Obj , RHS_Obj);
   ThomasAlgorithm::BackSubstitution(A3_Obj , RHS_Obj, Solution);
}
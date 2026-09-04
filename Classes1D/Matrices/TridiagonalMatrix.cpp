#include "Matrices/TridiagonalMatrix.hpp"

TridiagonalMatrix::TridiagonalMatrix(const Mesh1D& Mesh1D_Obj)
    :LowerDiagonalValues(Mesh1D_Obj.GetNumOfNodes(), 0.0),
     MainDiagonalValues(Mesh1D_Obj.GetNumOfNodes(), 0.0),
     UpperDiagonalValues(Mesh1D_Obj.GetNumOfNodes(), 0.0)
{
    
}

// IMPORTANT:
// LowerDiagonalValues[0] and UpperDiagonalValues[numNodes - 1] are unused
// placeholder elements. They only exist so that L[i], d[i], and u[i]
// can all be accessed using the same row index i.
//
// However, UpperDiagonalValues[0] and LowerDiagonalValues[numNodes - 1]
// are NOT unused:
//   - UpperDiagonalValues[0] is the coefficient of U[1] in the first row.
//   - LowerDiagonalValues[numNodes - 1] is the coefficient of U[numNodes - 2]
//     in the last row.
//
// For Dirichlet boundary conditions, the first and last matrix rows must
// directly enforce:
//   U[0] = leftBoundaryValue
//   U[numNodes - 1] = rightBoundaryValue
//
// Therefore, the boundary rows must be identity rows:
//   First row: d[0] = 1.0 and u[0] = 0.0
//   Last row:  L[numNodes - 1] = 0.0 and d[numNodes - 1] = 1.0
//
// The Crank-Nicolson coefficients L = -r/2, d = 1+r, and u = -r/2
// are valid only for the internal rows i = 1, ..., numNodes - 2.
// Do not assign -r/2 to u[0] or L[numNodes - 1], because doing so
// couples the prescribed boundary values to their neighboring nodes
// and prevents the matrix from directly enforcing the Dirichlet BCs.

std::size_t TridiagonalMatrix::SizeM() const
{
    return MainDiagonalValues.size();
}
std::size_t TridiagonalMatrix::SizeU()const
{
    return UpperDiagonalValues.size();
}

std::size_t TridiagonalMatrix::SizeL()const
{
    return LowerDiagonalValues.size();
}


double TridiagonalMatrix::GetMidValue(std::size_t Index_) const
{
    return MainDiagonalValues.at(Index_);
}
double TridiagonalMatrix::GetUValue(std::size_t Index_) const
{
    return UpperDiagonalValues.at(Index_);
}
double TridiagonalMatrix::GetLValue(std::size_t Index_) const
{
    return LowerDiagonalValues.at(Index_);
}



void TridiagonalMatrix::SetMValue(
    std::size_t Index_,
    double Value_)
{
    MainDiagonalValues.at(Index_) = Value_;
}

void TridiagonalMatrix::SetUValue(
    std::size_t Index_,
    double Value_)
{
    UpperDiagonalValues.at(Index_) = Value_;
}

void TridiagonalMatrix::SetLValue(
    std::size_t Index_,
    double Value_)
{
    LowerDiagonalValues.at(Index_) = Value_;
}






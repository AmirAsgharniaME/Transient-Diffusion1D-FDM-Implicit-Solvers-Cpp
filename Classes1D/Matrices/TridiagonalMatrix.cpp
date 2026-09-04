#include "Classes1D/Matrices/TridiagonalMatrix.hpp"

TridiagonalMatrix::TridiagonalMatrix(        
    SchemeName ActiveScheme,
    const Mesh1D& Mesh1D_Obj,
    SolutionParameters<double> DiffNumber_Obj,
    double Default_Initial_Value_ )
    :LowerDiagonalValues(Mesh1D_Obj.GetNumOfNodes(), Default_Initial_Value_),
     MainDiagonalValues(Mesh1D_Obj.GetNumOfNodes(), Default_Initial_Value_),
     UpperDiagonalValues(Mesh1D_Obj.GetNumOfNodes(), Default_Initial_Value_),
     r(DiffNumber_Obj.GetValue()),
     numNodes(Mesh1D_Obj.GetNumOfNodes())
{
    if (ActiveScheme == SchemeName::Laasonen)
    {
        // ================
        // Lassonen Method
        // ================
        // EQ : -r * U[i - 1][n + 1] + (1.0 + 2.0 * r) * U[i][n + 1] - r * U[i + 1][n + 1] = RHS(U_n)
        // RHS(U_n) = U[i][n]
        // Create : 
        // *LowerDiagonalValues For CrankNicolson Method
        // *MainDiagonalValues For CrankNicolson Method
        // *UpperDiagonalValues For CrankNicolson Method
        // ========================================

        const std::size_t numNodes = Mesh1D_Obj.GetNumOfNodes();

        // Iterate through the interior nodes
        for (std::size_t i = 1; i < numNodes - 1; ++i)
        {
            LowerDiagonalValues[i] = -r; //A[i][i - 1] = -r;

            MainDiagonalValues[i] = 1.0 + 2.0 * r; //A[i][i] = 1.0 + 2.0 * r;

            UpperDiagonalValues[i] = -r; //A[i][i + 1] = -r;
        }
    } 

    // ====================
    // CrankNicolson Method
    // ====================
    // EQ  : -(r / 2) * U[i - 1]^(n + 1) + (1 + r) * U[i]^(n + 1) - (r / 2) * U[i + 1]^(n + 1) = RHS(U_n)
    // RHS(U_n)=(r / 2) * U[i - 1]^n + (1 - r) * U[i]^n + (r / 2) * U[i + 1]^n
    // Create : 
    //     *LowerDiagonalValues For CrankNicolson Method
    //     *MainDiagonalValues For CrankNicolson Method
    //     *UpperDiagonalValues For CrankNicolson Method
    // =================================================

    if (ActiveScheme == SchemeName::CrankNicolson)
    {
        const std::size_t numNodes = Mesh1D_Obj.GetNumOfNodes();

        // Iterate through the interior nodes
        for (std::size_t i = 1 ; i < numNodes-1; ++i)
        {
            //LowerDiagonalValues For CrankNicolson Method
            LowerDiagonalValues[i] = -r / 2.0; // A[i][i - 1] = -r / 2.0;
            //MainDiagonalValues For CrankNicolson Method
            MainDiagonalValues[i] = 1.0 + r; //A[i][i] = 1.0 + r;
            //UpperDiagonalValues For CrankNicolson Method
            UpperDiagonalValues[i] = -r / 2.0; //A[i][i + 1] = -r / 2.0;
        }
    }
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

void TridiagonalMatrix::SetBoundaryConditions()
{
    //Apply Dirichlet boundary condition To MainDiagonal Values
    MainDiagonalValues[0] = 1.0;
    MainDiagonalValues[numNodes - 1] = 1.0;
   
    //Apply Dirichlet boundary condition To UpperDiagonal Values
    UpperDiagonalValues[0]= 0.0; //for Diriclet BC
    UpperDiagonalValues[ numNodes -1 ] = 0.0 ; // unused: just for Even index
    
    LowerDiagonalValues[0] = 0.0 ; // unused: just for Even index
    LowerDiagonalValues[numNodes -1] = 0.0 ; //for Diriclet BC
   
}




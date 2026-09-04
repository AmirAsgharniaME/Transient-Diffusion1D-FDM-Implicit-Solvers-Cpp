#include "Classes1D/Matrices/CoefficientMatrix.hpp"

CoefficientMatrix::CoefficientMatrix(
    Scheme ActiveScheme,
    const Mesh1D& Mesh1D_Obj,
    const SolutionParameters<double> DiffNumber_Obj)
    :nrows(Mesh1D_Obj.GetNumOfNodes()),
     ncols(Mesh1D_Obj.GetNumOfNodes()),
     AValues(nrows, std::vector<double>(ncols, 0.0)),
     r(DiffNumber_Obj.GetValue())
    {
    // ================
    // Lassonen Method
    // =================================================================================
    // EQ : -r * U[i - 1][n + 1] + (1.0 + 2.0 * r) * U[i][n + 1] - r * U[i + 1][n + 1] = RHS(U_n)
    // RHS(U_n) = U[i][n]
    // Create A for Laasonen Method
    // =================================================================================
        if (ActiveScheme == Scheme::Laasonen)
        {
            const std::size_t numNodes = Mesh1D_Obj.GetNumOfNodes();

            // Iterate through the interior nodes
            for (std::size_t i = 1; i < numNodes - 1; ++i)
            {
                AValues[i][i - 1] = -r; //A[i][i - 1] = -r;

                AValues[i][i] = 1.0 + 2.0 * r; //A[i][i] = 1.0 + 2.0 * r;

                AValues[i][i + 1] = -r; //A[i][i + 1] = -r;
            }
        }

    // =====================================================================================
    // CrankNicolson Method
    // =====================================================================================
    //  EQ : -(r / 2) * U[i - 1]^(n + 1) + (1 + r) * U[i]^(n + 1) - (r / 2) * U[i + 1]^(n + 1)
    //  = (r / 2) * U[i - 1]^n + (1 - r) * U[i]^n + (r / 2) * U[i + 1]^n
    //  Create A for CrankNicolson Method  
    // =====================================================================================
        if (ActiveScheme == Scheme::CrankNicolson)
        {
            const std::size_t numNodes = Mesh1D_Obj.GetNumOfNodes();

            // Iterate through the interior nodes
            for (std::size_t i = 1; i < numNodes - 1; ++i)
            {
                AValues[i][i - 1] = -r / 2.0; 

                AValues[i][i] = 1.0 + r; 

                AValues[i][i + 1] = -r / 2.0; 
            }
        }
        
    }


std::size_t CoefficientMatrix::GetNumRows() const
{
    return nrows;
}

std::size_t CoefficientMatrix::GetNumColumns() const
{
    return ncols;
}

double CoefficientMatrix::GetValue(
    std::size_t RowIndex,
    std::size_t ColumnIndex) const
{
    return AValues[RowIndex][ColumnIndex];
}

void CoefficientMatrix::SetValue(
    std::size_t RowIndex,
    std::size_t ColumnIndex,
    double Value_)
{
    AValues[RowIndex][ColumnIndex] = Value_;
}

void CoefficientMatrix::Print() const
{
    std::cout << std::fixed << std::setprecision(4);

    for (std::size_t j = 0; j < nrows; ++j)
    {
        for (std::size_t i = 0; i < ncols; ++i)
        {
            std::cout << std::setw(15) << AValues[j][i];
        }

        std::cout << '\n';
    }
}

void CoefficientMatrix::SetBoundaryConditions()
{
    
AValues[0][0] = 1.0;
AValues[nrows - 1][nrows - 1] = 1.0;

}
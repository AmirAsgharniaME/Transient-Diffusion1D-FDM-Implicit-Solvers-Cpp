#include "Matrices/CoefficientMatrix.hpp"

CoefficientMatrix::CoefficientMatrix(const Mesh1D& Mesh1D_Obj)
    :nrows(Mesh1D_Obj.GetNumOfNodes()),
     ncols(Mesh1D_Obj.GetNumOfNodes()),
     AValues(nrows, std::vector<double>(ncols, 0.0)) 
    {
        
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

void CoefficientMatrix::PrintMatrix() const
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


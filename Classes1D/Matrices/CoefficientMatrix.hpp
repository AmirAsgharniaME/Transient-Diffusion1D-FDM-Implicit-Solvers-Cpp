#ifndef CoefficientMatrix_HPP
#define CoefficientMatrix_HPP
#include <iostream>
#include <iomanip>
#include "Mesh1D/Mesh1D.hpp"
#include <string>
#include <vector>
#include <cstddef>

class CoefficientMatrix
{
public:
    CoefficientMatrix(const Mesh1D& Mesh1D_Obj);

[[nodiscard]] std::size_t GetNumRows() const;   //[[nodiscard]] : The output value of this function should not be left unused.
[[nodiscard]] std::size_t GetNumColumns() const;
[[nodiscard]] double GetValue(std::size_t RowIndex, std::size_t ColumnIndex) const;

void SetValue(
    std::size_t RowIndex,
    std::size_t ColumnIndex,
    double Value_);

 void PrintMatrix() const;


private:
    std::size_t nrows;
    std::size_t ncols;
    std::vector<std::vector<double>> AValues;
    
};

#endif

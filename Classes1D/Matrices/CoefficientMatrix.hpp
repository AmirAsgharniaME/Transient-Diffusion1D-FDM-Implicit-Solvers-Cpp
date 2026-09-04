#pragma once
#include <iostream>
#include <iomanip>
#include "Classes1D/Mesh1D/Mesh1D.hpp"
#include "Classes1D/SolverSettings/SolutionParameters.hpp"
#include <string>
#include <vector>
#include <cstddef>
#include "Classes1D/Boundaries/Boundary.hpp"

enum class  Scheme
{
    Laasonen,
    CrankNicolson
};

class CoefficientMatrix
{
public:
    CoefficientMatrix(
    Scheme ActiveScheme,
    const Mesh1D& Mesh1D_Obj,
    const SolutionParameters<double> DiffNumber_Obj);

[[nodiscard]] std::size_t GetNumRows() const;   //[[nodiscard]] : The output value of this function should not be left unused.
[[nodiscard]] std::size_t GetNumColumns() const;
[[nodiscard]] double GetValue(std::size_t RowIndex, std::size_t ColumnIndex) const;

void SetValue(
    std::size_t RowIndex,
    std::size_t ColumnIndex,
    double Value_);

 void Print() const;
 void SetBoundaryConditions();


private:
    std::size_t nrows;
    std::size_t ncols;
    std::vector<std::vector<double>> AValues;
     double r;
};


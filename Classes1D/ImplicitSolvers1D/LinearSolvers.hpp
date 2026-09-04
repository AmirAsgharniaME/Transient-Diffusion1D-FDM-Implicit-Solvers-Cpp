#ifndef LinearSolvers_HPP
#define LinearSolvers_HPP

#include <cstddef>
#include <iostream>
#include "Field1D/Field1D.hpp"
#include "Matrices/CoefficientMatrix.hpp"
#include "Matrices/TridiagonalMatrix.hpp"

#include "LinearSolvers/GaussianElimination.hpp"
#include "LinearSolvers/ThomasAlgorithm.hpp"

//#include  "CoefficientMatrix/TridiagonalMatrix.hpp"

//#include "LinearSolvers/ThomasAlgorithm.hpp"


class LinearSolvers
{
public:
    LinearSolvers() =delete;
    
   static void Solve_Gaussian_Elimination(CoefficientMatrix& A, RHS& RHS_Obj, Field1D& Solution);
   static void Solve_Thomas_Algorithm(TridiagonalMatrix& A3_Obj, RHS& RHS_Obj, Field1D& Solution);
};
#endif //ImplicitSolvers1D_HPP
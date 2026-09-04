#pragma once


#include "Classes1D/Field1D/Field1D.hpp"
#include  "Classes1D/Matrices/TridiagonalMatrix.hpp"
#include  "Classes1D/Matrices/RHS.hpp"
 
#include <cstddef>

class ThomasAlgorithm
{
public:
    // (Stateless Class)
    ThomasAlgorithm() = delete;


    static void Solve_nPlus1(
        TridiagonalMatrix& TridiagonalMatrix,
        RHS& RHS_Obj,
        Field1D& Solution);

//private:
    // (Forward Sweep)
    static void ForwardElimination(
        TridiagonalMatrix& TridiagonalMatrix,
        RHS& RHS_Obj);

    // (Back Substitution)
    static void BackSubstitution(
        const TridiagonalMatrix& TridiagonalMatrix,
        const RHS& RHS_Obj,
        Field1D& Solution);
};

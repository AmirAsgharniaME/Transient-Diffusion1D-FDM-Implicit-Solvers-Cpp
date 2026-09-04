#ifndef ThomasAlgorithm_HPP
#define ThomasAlgorithm_HPP

#include "Field1D/Field1D.hpp"
#include  "Matrices/TridiagonalMatrix.hpp"
#include  "Matrices/RHS.hpp"
 
#include <cstddef>

class ThomasAlgorithm
{
public:
    // (Stateless Class)
    ThomasAlgorithm() = delete;


    static void Solve(
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

#endif

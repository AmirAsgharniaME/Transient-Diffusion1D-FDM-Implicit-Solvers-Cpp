#ifndef GaussianElimination_HPP
#define GaussianElimination_HPP

#include "Field1D/Field1D.hpp"
#include "Matrices/CoefficientMatrix.hpp"
#include "Matrices/RHS.hpp"

#include <cstddef>
//Gaussian Elimination with Partial Row Pivoting followed by Back Substitution
class GaussianElimination
{
public:

 GaussianElimination() = delete;   //Creating an object of this class is prohibited.
    //A Solution = RHS
    static void Solve(
        CoefficientMatrix& A,
        RHS& RHS,
        Field1D& Solution);

    static void ForwardElimination(
        CoefficientMatrix& A,
        RHS& RHS);
/*
 * Back Substitution Formula:
 * 
 *   x[i] = ( b[i] - Sum_{j = i+1}^{n-1} ( a[i][j] * x[j] ) ) / a[i][i]
 * 
 * Execution Steps:
 *   1. Initialize 'Sum' as b[i] (RHS value).
 *   2. Subtract the product of known variables: Sum -= a[i][j] * x[j] (for j from i+1 to n-1).
 *   3. Calculate the final value: x[i] = Sum / a[i][i].
 */


    static void BackSubstitution(
        const CoefficientMatrix& A,
        const RHS& RHS,
        Field1D& Solution);

    [[nodiscard]] static std::size_t FindPivotRow(
        const CoefficientMatrix& A,
        std::size_t PivotIndex);

    static void SwapRows(
        CoefficientMatrix& A,
        RHS& RHS,
        std::size_t FirstRow,
        std::size_t SecondRow);
};

#endif

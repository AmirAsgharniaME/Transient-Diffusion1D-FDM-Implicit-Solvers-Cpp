#include "LinearSolvers/GaussianElimination.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

void GaussianElimination::Solve(
    CoefficientMatrix& A,
    RHS& RHS,
    Field1D& Solution)
{
    const std::size_t NumRows = A.GetNumRows();
    const std::size_t NumCols = A.GetNumColumns();

    if (NumRows != NumCols)
    {
        throw std::invalid_argument("Solve: matrix A must be square.");
    }

    if (RHS.GetSize() != NumRows || Solution.Size() != NumRows)
    {
        throw std::invalid_argument("Solve: size of RHS and Solution must match matrix size.");
    }

    ForwardElimination(A, RHS);
    BackSubstitution(A, RHS, Solution);
}

void GaussianElimination::ForwardElimination(
    CoefficientMatrix& A,
    RHS& RHS)
{
    const std::size_t NumRows = A.GetNumRows();
    const std::size_t NumCols = A.GetNumColumns();

    if (NumRows != NumCols)
    {
        throw std::invalid_argument("ForwardElimination: matrix A must be square.");
    }

    if (RHS.GetSize() != NumRows)
    {
        throw std::invalid_argument("ForwardElimination: size of RHS must match matrix size.");
    }

    if (NumRows == 0)
    {
        return;
    }

    for (std::size_t PivotIndex = 0; PivotIndex < NumRows - 1; ++PivotIndex)
    {
        const std::size_t PivotRow = FindPivotRow(A, PivotIndex);
        const double PivotValue = A.GetValue(PivotRow, PivotIndex);

        if (std::fabs(PivotValue) < 1.0e-14)
        {
            throw std::runtime_error("ForwardElimination: singular or nearly singular matrix detected.");
        }

        if (PivotRow != PivotIndex)
        {
            SwapRows(A, RHS, PivotIndex, PivotRow);
        }

        const double DiagonalValue = A.GetValue(PivotIndex, PivotIndex);

        for (std::size_t RowIndex = PivotIndex + 1; RowIndex < NumRows; ++RowIndex)
        {
            const double Factor = A.GetValue(RowIndex, PivotIndex) / DiagonalValue;

            A.SetValue(RowIndex, PivotIndex, 0.0);

            for (std::size_t ColIndex = PivotIndex + 1; ColIndex < NumCols; ++ColIndex)
            {
                const double NewValue =
                    A.GetValue(RowIndex, ColIndex) - Factor * A.GetValue(PivotIndex, ColIndex);
                A.SetValue(RowIndex, ColIndex, NewValue);
            }

            const double NewRHS =
                RHS.GetValue(RowIndex) - Factor * RHS.GetValue(PivotIndex);
            RHS.SetValue(RowIndex, NewRHS);
        }
    }

    if (std::fabs(A.GetValue(NumRows - 1, NumRows - 1)) < 1.0e-14)
    {
        throw std::runtime_error("ForwardElimination: zero pivot found on last row.");
    }
}

void GaussianElimination::BackSubstitution(
    const CoefficientMatrix& A,
    const RHS& RHS,
    Field1D& Solution)
{
    const std::size_t NumRows = A.GetNumRows();
    const std::size_t NumCols = A.GetNumColumns();

    if (NumRows != NumCols)
    {
        throw std::invalid_argument("BackSubstitution: matrix A must be square.");
    }

    if (RHS.GetSize() != NumRows || Solution.Size() != NumRows)
    {
        throw std::invalid_argument("BackSubstitution: size of RHS and Solution must match matrix size.");
    }
    //Standard and safe way to iterate backward using an unsigned index.
    // Safely iterates backward without causing unsigned index underflow.
    for (std::size_t RowIndex = NumRows; RowIndex-- > 0;)
    {
        double Sum = RHS.GetValue(RowIndex);

        // For the last row (RowIndex = NumRows - 1), this loop is skipped because ColIndex starts at NumCols.
        // Skip summation for the last row since there are no post-diagonal elements to substitute.
        for (std::size_t ColIndex = RowIndex + 1; ColIndex < NumCols; ++ColIndex)
        {
            Sum -= A.GetValue(RowIndex, ColIndex) * Solution.GetValue(ColIndex);
        }

        const double DiagonalValue = A.GetValue(RowIndex, RowIndex);

        if (std::fabs(DiagonalValue) < 1.0e-14)
        {
            throw std::runtime_error("BackSubstitution: zero diagonal entry detected.");
        }

        Solution.SetValue(RowIndex, Sum / DiagonalValue);
    }
}

std::size_t GaussianElimination::FindPivotRow(
    const CoefficientMatrix& A,
    std::size_t PivotIndex)
{
    const std::size_t NumRows = A.GetNumRows();

    if (PivotIndex >= NumRows)
    {
        throw std::out_of_range("FindPivotRow: PivotIndex is out of range.");
    }

    std::size_t PivotRow = PivotIndex;
    double MaxAbsValue = std::fabs(A.GetValue(PivotIndex, PivotIndex));

    for (std::size_t RowIndex = PivotIndex + 1; RowIndex < NumRows; ++RowIndex)
    {
        const double CurrentAbsValue = std::fabs(A.GetValue(RowIndex, PivotIndex));
        if (CurrentAbsValue > MaxAbsValue)
        {
            MaxAbsValue = CurrentAbsValue;
            PivotRow = RowIndex;
        }
    }

    return PivotRow;
}

void GaussianElimination::SwapRows(
    CoefficientMatrix& A,
    RHS& RHS,
    std::size_t FirstRow,
    std::size_t SecondRow)
{
    const std::size_t NumCols = A.GetNumColumns();

    if (FirstRow >= A.GetNumRows() || SecondRow >= A.GetNumRows())
    {
        throw std::out_of_range("SwapRows: row index is out of range.");
    }

    if (RHS.GetSize() != A.GetNumRows())
    {
        throw std::invalid_argument("SwapRows: size of RHS must match matrix size.");
    }

    if (FirstRow == SecondRow)
    {
        return;
    }

    for (std::size_t ColIndex = 0; ColIndex < NumCols; ++ColIndex)
    {
        const double Temp = A.GetValue(FirstRow, ColIndex);
        A.SetValue(FirstRow, ColIndex, A.GetValue(SecondRow, ColIndex));
        A.SetValue(SecondRow, ColIndex, Temp);
    }

    const double TempRHS = RHS.GetValue(FirstRow);
    RHS.SetValue(FirstRow, RHS.GetValue(SecondRow));
    RHS.SetValue(SecondRow, TempRHS);
}

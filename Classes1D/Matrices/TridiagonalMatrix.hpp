#pragma once


#include "Classes1D/Mesh1D/Mesh1D.hpp"
#include <string>
#include <vector>
#include <cstddef>
#include "Classes1D/SolverSettings/SolutionParameters.hpp"

enum class  SchemeName
{
    Laasonen,
    CrankNicolson
};

class TridiagonalMatrix
{
public:
    TridiagonalMatrix(
        SchemeName ActiveScheme,
        const Mesh1D& Mesh1D_Obj,
        SolutionParameters<double> DiffNumber_Obj,
        double Default_Initial_Value_ = 0.0);





std::size_t SizeM() const;
std::size_t SizeU() const;
std::size_t SizeL() const;

double GetMidValue(std::size_t Index_) const;
double GetUValue(std::size_t Index_) const;
double GetLValue(std::size_t Index_) const;


void SetMValue(
    std::size_t Index_,
    double Value_);

void SetUValue(
    std::size_t Index_,
    double Value_);

void SetLValue(
    std::size_t Index_,
    double Value_);

void SetBoundaryConditions();


private:
    std::vector<double> LowerDiagonalValues;
    std::vector<double> MainDiagonalValues;
    std::vector<double> UpperDiagonalValues;
    double r;
    const std::size_t numNodes;
};



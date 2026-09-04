#ifndef TridiagonalMatrix_HPP
#define TridiagonalMatrix_HPP

#include "Mesh1D/Mesh1D.hpp"
#include <string>
#include <vector>
#include <cstddef>

class TridiagonalMatrix
{
public:
    TridiagonalMatrix(const Mesh1D& Mesh1D_Obj);





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


private:
    std::vector<double> LowerDiagonalValues;
    std::vector<double> MainDiagonalValues;
    std::vector<double> UpperDiagonalValues;
};

#endif //TridiagonalCoefficientMatrix_HPP

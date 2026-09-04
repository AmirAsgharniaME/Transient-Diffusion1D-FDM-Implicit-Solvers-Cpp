#ifndef RHS_HPP
#define RHS_HPP
#include <iostream>
#include <iomanip>
#include <cstddef>
#include <vector>
#include <string>
#include "Field1D/Field1D.hpp"
#include "SolverSettings/SolutionParameters.hpp"

class RHS
{


private:
    std::vector<double> RHSValues;
    double DiffNumer;
  


public:

    //RHS For Lassonen Method
    explicit RHS(
        const Field1D& Field1D_Obj,
        double Default_Initial_Value_ = 0.0
    );


    //RHS For CrankNicolson Method
    explicit RHS(
        const Field1D& Field1D_Obj,
        const SolutionParameters<double>& DiffNumer_Obj,
        double Default_Initial_Value_ = 0.0      
    );

    
    /*****Public Setter*****/ 
    void SetValue(std::size_t PIndex, double Value_);
    void SetAllSingleValue(double Value_);

    /*****Public Getters*****/ 
     std::size_t GetSize() const;
     double GetValue(std::size_t PIndex) const;


    
    



};



#endif //RHS_HPP
